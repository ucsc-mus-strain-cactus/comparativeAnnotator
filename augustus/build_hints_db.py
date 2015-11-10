"""
This program takes a series of BAM files and converts them to a Augustus hints database.
"""

import sys
import os
import pysam
import argparse
from pyfaidx import Fasta
from lib.general_lib import format_ratio, get_tmp
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


def main_hints_fn(target, bam_paths, db_path, genome, genome_fasta):
    filtered_bam_tree = TempFileTree(get_tmp(target, global_dir=True, name="filter_file_tree"))
    for bam_path in bam_paths:
        paired = "--paired --pairwiseAlignments" if bam_is_paired(bam_path) is True else ""
        sam_handle = pysam.Samfile(bam_path)
        for reference in sam_handle.references:
            out_filter = filtered_bam_tree.getTempFile(suffix=".bam")
            target.addChildTargetFn(sort_by_name, memory=8 * 1024 ** 3, cpu=2, 
                                    args=[bam_path, reference, out_filter, paired])
    target.setFollowOnTargetFn(build_hints, args=[filtered_bam_tree, genome, db_path, genome_fasta])


def sort_by_name(target, bam_path, references, out_filter, paired):
    assert references 
    references = " ".join(references)
    tmp_filtered = get_tmp(target, name=".filtered.bam")
    cmd = "samtools view -b {} {} | samtools sort -O bam -T {} -n - | filterBam --uniq {} --in /dev/stdin --out {}"
    cmd = cmd.format(bam_path, references, get_tmp(target), paired, tmp_filtered)
    system(cmd)
    cmd2 = "samtools sort -O bam -T {} {} > {}".format(get_tmp(target), tmp_filtered, out_filter)
    system(cmd2)
    system("samtools index {}".format(out_filter))


def build_hints(target, filtered_bam_tree, genome, db_path, genome_fasta):
    bam_files = [x for x in file_tree.listFiles() if x.endswith("bam")]
    intron_hints_tree = TempFileTree(get_tmp(target, global_dir=True, name="intron_hints_tree"))
    exon_hints_tree = TempFileTree(get_tmp(target, global_dir=True, name="exon_hints_tree"))
    for bam_file in bam_files:
        intron_hints_path = intron_hints_tree.getTempFile(suffix=".intron.gff")
        target.addChildTargetFn(build_intron_hints, memory=8 * 1024 ** 3, cpu=2, args=[bam_file, intron_hints_path])
        exon_hints_path = intron_hints_tree.getTempFile(suffix=".exon.gff")
        target.addChildTargetFn(build_exon_hints, memory=8 * 1024 ** 3, cpu=2, args=[bam_file, exon_hints_path])
    target.setFollowOnTargetFn(cat_hints, args[intron_hints_tree, exon_hints_tree, genome, db_path, genome_fasta])


def build_exon_hints(target, bam_file, exon_gff_path):
    cmd = ('bam2wig {} | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep '
           '--UCSC=/dev/null --radius=4.5 --pri=4 --strand="." > {}')
    cmd = cmd.format(bam_file, exon_gff_path)
    system(cmd)


def build_intron_hints(target, bam_file, intron_hints_path):
    cmd = "bam2hints --intronsonly --in {} --out {}"
    cmd = cmd.format(bam_file, intron_hints_path)
    system(cmd)


def cat_hints(target, intron_hints_tree, exon_hints_tree, genome, db_path, genome_fasta):
    all_gffs = intron_hints_tree.listFiles() + exon_hints_tree.listFiles()
    hints = get_tmp(target, global_dir=True, name="combined_sorted_hints.gff")
    cmd = "cat {} | sort -n -k4,4 | sort -s -n -k5,5 | sort -s -n -k3,3 | sort -s -k1,1 | join_mult_hints.pl > {}"
    cmd = cmd.format(all_gffs, hints)
    system(cmd)
    target.setFollowOnTargetFn(load_db, args=[hints, db_path, genome, genome_fasta])


def load_db(target, hints, db_path, genome, genome_fasta):
    cmd = "load2sqlitedb --noIdx --species={} --dbaccess={} {}"
    fa_cmd = cmd.format(genome, db_path, genome_fasta)
    system(fa_cmd)
    hints_cmd = cmd.format(genome, db_path, hints)
    system(hints_cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--filterTissues", nargs="+")
    parser.add_argument("--filterCenters", nargs="+")
    bamfiles = parser.add_mutually_exclusive_group(required=True)
    bamfiles.add_argument("--bamFiles", nargs="+", help="bamfiles being used", dest="bams")
    bamfiles.add_argument("--bamFofn", help="File containing list of bamfiles", dest="bams")
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
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
                                  args=[args.bams, args.database, args.genome, args.fasta]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from augustus.build_hints_db import *
    main()
