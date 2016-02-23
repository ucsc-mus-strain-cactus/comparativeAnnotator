"""
This program takes a series of BAM files and converts them to a Augustus hints database.
"""

import sys
import time
import os
import pysam
import argparse
import itertools
import subprocess
from pyfaidx import Fasta
from lib.general_lib import format_ratio, get_tmp, mkdir_p
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, popenCatch, getRandomAlphaNumericString, catFiles, TempFileTree


def bam_is_paired(path, num_reads=100000, paired_cutoff=0.75):
    """
    Infers the paired-ness of a bam file.
    """
    sam = pysam.Samfile(path)
    count = 0
    for i, rec in enumerate(sam):
        if rec.is_paired:
            count += 1
        if i == num_reads:
            break
    if format_ratio(count, num_reads) > 0.75:
        return True
    elif format_ratio(count, num_reads) < 1 - paired_cutoff:
        return False
    else:
        raise RuntimeError("Unable to infer pairing from bamfile {}".format(path))


def group_references(sam_handle, num_bases=20 ** 7, max_seqs=100):
    """
    Group up references by num_bases, unless that exceeds max_seqs.
    """
    name_iter = itertools.izip(*[sam_handle.references, sam_handle.lengths])
    name, size = name_iter.next()
    this_bin = [name]
    bin_base_count = size
    num_seqs = 1
    for name, size in name_iter:
        bin_base_count += size
        num_seqs += 1
        if bin_base_count >= num_bases or num_seqs > max_seqs:
            yield this_bin
            this_bin = [name]
            bin_base_count = size
            num_seqs = 1
        else:
            this_bin.append(name)


def main_hints_fn(target, bam_paths, db_path, genome, genome_fasta, hints_dir):
    """
    Main driver function. Loops over each BAM, inferring paired-ness, then passing each BAM with one chromosome name
    for filtering. Each BAM will remain separated until the final concatenation and sorting of the hint gffs.
    """
    filtered_bam_tree = TempFileTree(get_tmp(target, global_dir=True, name="filter_file_tree"))
    for bam_path in bam_paths:
        paired = "--paired --pairwiseAlignments" if bam_is_paired(bam_path) is True else ""
        sam_handle = pysam.Samfile(bam_path)
        for references in group_references(sam_handle):
            out_filter = filtered_bam_tree.getTempFile(suffix=".bam")
            target.addChildTargetFn(sort_by_name, memory=8 * 1024 ** 3, cpu=2, 
                                    args=[bam_path, references, out_filter, paired])
    target.setFollowOnTargetFn(build_hints, args=[filtered_bam_tree, genome, db_path, genome_fasta, hints_dir])


def sort_by_name(target, bam_path, references, out_filter, paired):
    """
    Slices out a chromosome from a BAM, re-sorts by name, filters the reads, then re-sorts by position.
    filterBam does weird things when piped to stdout, so I don't do that.
    """
    references = " ".join(references)
    tmp_filtered = get_tmp(target, name=".filtered.bam")
    cmd = "samtools view -b {} {} | samtools sort -O bam -T {} -n - | filterBam --uniq {} --in /dev/stdin --out {}"
    cmd = cmd.format(bam_path, references, get_tmp(target), paired, tmp_filtered)
    system(cmd)
    cmd2 = "samtools sort -O bam -T {} {} > {}".format(get_tmp(target), tmp_filtered, out_filter)
    system(cmd2)
    system("samtools index {}".format(out_filter))


def build_hints(target, filtered_bam_tree, genome, db_path, genome_fasta, hints_dir):
    """
    Driver function for hint building. Builts intron and exon hints, then calls cat_hints to do final concatenation
    and sorting.
    """
    bam_files = [x for x in filtered_bam_tree.listFiles() if x.endswith("bam")]
    intron_hints_tree = TempFileTree(get_tmp(target, global_dir=True, name="intron_hints_tree"))
    exon_hints_tree = TempFileTree(get_tmp(target, global_dir=True, name="exon_hints_tree"))
    for bam_file in bam_files:
        intron_hints_path = intron_hints_tree.getTempFile(suffix=".intron.gff")
        target.addChildTargetFn(build_intron_hints, memory=8 * 1024 ** 3, cpu=2, args=[bam_file, intron_hints_path])
        exon_hints_path = exon_hints_tree.getTempFile(suffix=".exon.gff")
        target.addChildTargetFn(build_exon_hints, memory=8 * 1024 ** 3, cpu=2, args=[bam_file, exon_hints_path])
    target.setFollowOnTargetFn(cat_hints, args=[intron_hints_tree, exon_hints_tree, genome, db_path, genome_fasta,
                                                hints_dir])


def build_exon_hints(target, bam_file, exon_gff_path):
    """
    Exon hints are built from peaks of coverage, using bam2wig and wig2hints.pl in the Augustus repo.
    """
    cmd = ('bam2wig {} | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep '
           '--UCSC=/dev/null --radius=4.5 --pri=4 --strand="." > {}')
    cmd = cmd.format(bam_file, exon_gff_path)
    system(cmd)


def build_intron_hints(target, bam_file, intron_hints_path):
    """
    Intron hints look for splice junctions in BAM files
    """
    cmd = "bam2hints --intronsonly --in {} --out {}"
    cmd = cmd.format(bam_file, intron_hints_path)
    system(cmd)


def cat_hints(target, intron_hints_tree, exon_hints_tree, genome, db_path, genome_fasta, hints_dir):
    """
    All intron and exon hint gff files are concatenated and then sorted.
    """
    all_gffs = intron_hints_tree.listFiles() + exon_hints_tree.listFiles()
    # we use a fofn to side-step cat's inability to handle long input strings
    gff_fofn = get_tmp(target, name="gff_fofn")
    with open(gff_fofn, "w") as outf:
        for x in all_gffs:
            outf.write(x + "\n")
    concat_hints = get_tmp(target, name="concat_hints")
    cat_cmd = "cat {} | xargs -n 50 cat >> {}".format(gff_fofn, concat_hints)
    system(cat_cmd)
    hints = os.path.join(hints_dir, genome + ".reduced_hints.gff")
    # TODO: this takes forever. Surely this can be merged into one better sort command
    cmd = "cat {} | sort -n -k4,4 | sort -s -n -k5,5 | sort -s -n -k3,3 | sort -s -k1,1 | join_mult_hints.pl > {}"
    cmd = cmd.format(concat_hints, hints)
    system(cmd)
    target.setFollowOnTargetFn(load_db, args=[hints, db_path, genome, genome_fasta])


def load_db(target, hints, db_path, genome, genome_fasta, timeout=30000, intervals=120):
    """
    Final database loading.
    NOTE: Once done on all genomes, you want to run load2sqlitedb --makeIdx --dbaccess ${db}
    """
    cmd = "load2sqlitedb --noIdx --species={} --dbaccess={} {}"
    fa_cmd = cmd.format(genome, db_path, genome_fasta)
    hints_cmd = cmd.format(genome, db_path, hints)
    def handle_concurrency(cmd, timeout, intervals, start_time=None):
        if start_time is None:
            start_time = time.time()
        elif time.time() - start_time >= timeout:
            raise RuntimeError("hints database still locked after {} seconds".format(timeout))
        p = subprocess.Popen(cmd, shell=True, bufsize=-1, stderr=subprocess.PIPE)
        _, ret = p.communicate()
        if p.returncode == 0:
            return 1
        elif p.returncode == 1 and "locked" in ret:
            time.sleep(intervals)
            handle_concurrency(cmd, timeout, intervals, start_time)
        else:
            raise RuntimeError(ret)
    mkdir_p(os.path.dirname(db_path))
    for cmd in [fa_cmd, hints_cmd]:
        ret = handle_concurrency(cmd, timeout, intervals)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--hintsDir", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--filterTissues", nargs="+")
    parser.add_argument("--filterCenters", nargs="+")
    bamfiles = parser.add_mutually_exclusive_group(required=True)
    bamfiles.add_argument("--bamFiles", nargs="+", help="bamfiles being used", dest="bams")
    bamfiles.add_argument("--bamFofn", help="File containing list of bamfiles", dest="bams")
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    # below is an ugly hack to filter out tissues/centers by checking for the words in the file path
    if not isinstance(args.bams, list):
        if not os.path.exists(args.bams):
            raise RuntimeError("ERROR: bamFofn does not exist.")
        bams = {x.rstrip() for x in open(args.bams)}
        args.bams = bams
    else:
        args.bams = set(args.bams)
    for to_remove_list in [args.filterTissues, args.filterCenters]:
        if isinstance(to_remove_list, list):
            to_remove = set()
            for x in to_remove_list:
                for b in args.bams:
                    if x in b:
                        to_remove.add(b)
            args.bams -= to_remove
    s = Stack(Target.makeTargetFn(main_hints_fn, memory=8 * 1024 ** 3,
                                  args=[args.bams, args.database, args.genome, args.fasta, args.hintsDir]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from augustus.build_hints_db import *
    main()
