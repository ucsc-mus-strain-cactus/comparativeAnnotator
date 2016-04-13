"""
This program takes a series of BAM files and converts them to a Augustus hints database.
"""
import sys
import os
import pysam
import time
import argparse
import itertools
import subprocess
from collections import defaultdict
from pyfasta import Fasta
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from pycbio.sys.mathOps import format_ratio
from pycbio.sys.fileOps import tmpFileGet, ensureDir, atomicInstall
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, popenCatch, getRandomAlphaNumericString, catFiles, TempFileTree


def get_tmp(target, global_dir=False, name=None):
    """
    Wrapper functionality for getting temp dirs from jobTree.
    """
    prefix = getRandomAlphaNumericString(10)
    name = "".join([prefix, name]) if name is not None else prefix
    if global_dir is False:
        return os.path.join(target.getLocalTempDir(), name)
    else:
        return os.path.join(target.getGlobalTempDir(), name)


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


def group_references(sam_handle, num_bases=10 ** 6, max_seqs=100):
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


def main_hints_fn(target, bam_paths, db, genome, genome_fasta, hints_file):
    """
    Main driver function. Loops over each BAM, inferring paired-ness, then passing each BAM with one chromosome name
    for filtering. Maintains a chromosome based map of bams so that they can be merged prior to hint construction.
    """
    filtered_bam_tree = TempFileTree(get_tmp(target, global_dir=True, name="filter_file_tree"))
    chrom_file_map = defaultdict(list)
    sam_handle = pysam.Samfile(bam_paths[0])
    references = list(group_references(sam_handle))
    for bam_path in bam_paths:
        paired = "--paired --pairwiseAlignments" if bam_is_paired(bam_path) is True else ""
        for group in references:
            out_filter = filtered_bam_tree.getTempFile(suffix=".bam")
            chrom_file_map[tuple(group)].append(out_filter)
            target.addChildTargetFn(sort_by_name, args=[bam_path, group, out_filter, paired])
    target.setFollowOnTargetFn(merge_bam_wrapper, args=[chrom_file_map, genome, db, genome_fasta, hints_file])


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


def merge_bam_wrapper(target, chrom_file_map, genome, db, genome_fasta, hints_file):
    """
    Wrapper to parallelize merging of bams by chromosome.
    """
    merged_bams = []
    for chrom, bam_files in chrom_file_map.iteritems():
        merged_path = get_tmp(target, global_dir=True, name='.combined.bam')
        target.addChildTargetFn(merge_bams, args=(bam_files, merged_path))
        merged_bams.append(merged_path)
    target.setFollowOnTargetFn(build_hints, args=[merged_bams, genome, db, genome_fasta, hints_file])


def merge_bams(target, bam_files, merged_path):
    fofn = get_tmp(target)
    with open(fofn, 'w') as outf:
        for x in bam_files:
            outf.write(x + "\n")
    cmd = 'samtools merge -b {} {}'.format(fofn, merged_path)
    system(cmd)


def build_hints(target, merged_bams, genome, db, genome_fasta, hints_file):
    """
    Driver function for hint building. Builts intron and exon hints, then calls cat_hints to do final concatenation
    and sorting.
    """
    if hints_file is None:
        hints_file = target.getGlobalTempDir()
    intron_hints_tree = TempFileTree(get_tmp(target, global_dir=True, name="intron_hints_tree"))
    exon_hints_tree = TempFileTree(get_tmp(target, global_dir=True, name="exon_hints_tree"))
    for merged_bam in merged_bams:
        intron_hints_path = intron_hints_tree.getTempFile(suffix=".intron.gff")
        target.addChildTargetFn(build_intron_hints, memory=8 * 1024 ** 3, args=[merged_bam, intron_hints_path])
        exon_hints_path = exon_hints_tree.getTempFile(suffix=".exon.gff")
        target.addChildTargetFn(build_exon_hints, memory=8 * 1024 ** 3, args=[merged_bam, exon_hints_path])
    target.setFollowOnTargetFn(cat_hints, args=[intron_hints_tree, exon_hints_tree, genome, db, genome_fasta,
                                                hints_file])


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


def cat_hints(target, intron_hints_tree, exon_hints_tree, genome, db, genome_fasta, hints_file):
    """
    All intron and exon hint gff files are concatenated and then sorted.
    """
    if os.path.dirname(hints_file) != '':
        ensureDir(os.path.dirname(hints_file))
    all_gffs = intron_hints_tree.listFiles() + exon_hints_tree.listFiles()
    # we use a fofn to side-step cat's inability to handle long input strings
    gff_fofn = get_tmp(target, name="gff_fofn")
    with open(gff_fofn, "w") as outf:
        for x in all_gffs:
            outf.write(x + "\n")
    concat_hints = get_tmp(target, name="concat_hints")
    cat_cmd = "cat {} | xargs -n 50 cat >> {}".format(gff_fofn, concat_hints)
    system(cat_cmd)
    hints_tmp = tmpFileGet()
    # TODO: this takes forever. Surely this can be merged into one better sort command
    cmd = "cat {} | sort -n -k4,4 | sort -s -n -k5,5 | sort -s -n -k3,3 | sort -s -k1,1 | join_mult_hints.pl > {}"
    cmd = cmd.format(concat_hints, hints_tmp)
    system(cmd)
    atomicInstall(hints_tmp, hints_file)
    target.setFollowOnTargetFn(load_db, args=[hints_file, db, genome, genome_fasta])


def load_db(target, hints_file, db, genome, genome_fasta, timeout=6000, intervals=30):
    """
    Final database loading.
    NOTE: Once done on all genomes, you want to run load2sqlitedb --makeIdx --dbaccess ${db}
    """
    cmd = "load2sqlitedb --noIdx --species={} --dbaccess={} {}"
    fa_cmd = cmd.format(genome, db, genome_fasta)
    hints_cmd = cmd.format(genome, db, hints_file)
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
    if os.path.dirname(db) != '':
        ensureDir(os.path.dirname(db))
    for cmd in [fa_cmd, hints_cmd]:
        ret = handle_concurrency(cmd, timeout, intervals)


def external_main(args):
    assert os.path.exists(args.fasta)
    assert all([os.path.exists(x) for x in args.bams])
    args.defaultMemory = 8 * 1024 ** 3
    s = Stack(Target.makeTargetFn(main_hints_fn, args=[args.bams, args.database, args.genome, args.fasta,
                                                       args.hintsFile]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--hintsFile", default=None)
    parser.add_argument("--fasta", required=True)
    bamfiles = parser.add_mutually_exclusive_group(required=True)
    bamfiles.add_argument("--bamFiles", nargs="+", help="bamfiles being used", dest="bams")
    bamfiles.add_argument("--bamFofn", help="File containing list of bamfiles", dest="bams")
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    # below is an ugly hack to filter out tissues/centers by checking for the words in the file path
    assert os.path.exists(args.fasta)
    if not isinstance(args.bams, list):
        if not os.path.exists(args.bams):
            raise RuntimeError("ERROR: bamFofn does not exist.")
        bams = {x.rstrip() for x in open(args.bams)}
        args.bams = bams
    else:
        args.bams = set(args.bams)
    args.defaultMemory = 8 * 1024 ** 3
    s = Stack(Target.makeTargetFn(main_hints_fn, args=[args.bams, args.database, args.genome, args.fasta,
                                                       args.hintsFile]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from comparativeAnnotator.augustus.build_hints_db import *
    main()
