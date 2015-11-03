"""
This program takes a series of BAM files and converts them to a Augustus hints database.
"""

import sys
import os
import pysam
import argparse
from pyfaidx import Fasta
from lib.general_lib import format_ratio
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


def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


def get_tmp(target, global_dir=False, name=None):
    name = getRandomAlphaNumericString(10) if name is None else name
    if global_dir is False:
        return os.path.join(target.getLocalTempDir(), name)
    else:
        return os.path.join(target.getGlobalTempDir(), name)


def filter_bam(target, path, file_tree):
    paired = "--paired --pairwiseAlignments" if bam_is_paired(path) is True else ""
    cmd = ("samtools sort -O bam -T {} -n {} | filterBam --uniq {} --in /dev/stdin --out /dev/stdout | samtools sort "
           "-O bam -T {} - > {}")
    bam_path = file_tree.getTempFile(suffix="bam")
    cmd = cmd.format(get_tmp(target), path, paired, get_tmp(target), bam_path)
    system(cmd)
    system("samtools index {}".format(bam_path))


def merge_bam_to_wig(target, file_tree, wig_path, db_path, genome, genome_fasta):
    bam_files = [x for x in file_tree.listFiles() if x.endswith("bam")]
    fofn_path = os.path.join(target.getGlobalTempDir(), "fofn")
    with open(fofn_path, "w") as fofn:
        for x in bam_files:
            fofn.write(x + "\n")
    cmd = "bamtools merge -list {} | bam2wig /dev/stdin > {}".format(fofn_path, wig_path)
    system(cmd)
    exon_gff_path = get_tmp(target, global_dir=True, name="exon_hints.gff")
    intron_gff_path = get_tmp(target, global_dir=True, name="intron_hints.gff")
    target.addChildTargetFn(wig_2_hints, memory=8 * 1024 ** 3, cpu=2,
                            args=[wig_path, exon_gff_path])
    target.addChildTargetFn(bam_2_hints_wrapper, args=[bam_files, intron_gff_path, genome_fasta])
    target.setFollowOnTargetFn(cat_hints, args=[exon_gff_path, intron_gff_path, db_path, genome, genome_fasta])


def wig_2_hints(target, wig_path, exon_gff_path):
    cmd = ('cat {} | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep '
           '--UCSC=/dev/null --radius=4.5 --pri=4 --strand="." > {}')
    cmd = cmd.format(wig_path, exon_gff_path)
    system(cmd)


def bam_2_hints_wrapper(target, bam_files, intron_gff_path, genome_fasta):
    intron_hints_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), "intron_hints_tree"))
    tmp_bam_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), "tmp_bam_tree"))
    fasta = Fasta(genome_fasta)
    if len(fasta.keys()) > 1000:
        keys = chunker(fasta.keys(), 250)
    else:
        keys = fasta.keys()
    for f in bam_files:
        for chroms in keys:
            target.addChildTargetFn(bam_2_hints, memory=8 * 1024 ** 3, cpu=2, 
                                    args=[f, intron_hints_tree, chroms])
    target.setFollowOnTargetFn(cat_bam_2_hints, memory=8 * 1024 ** 3, cpu=4,
                                args=[intron_hints_tree, intron_gff_path])


def bam_2_hints(target, bam_file, intron_hints_tree, chroms):
    chroms = " ".join(chroms)
    f = intron_hints_tree.getTempFile(suffix="hints")
    tmp = get_tmp(target, name="tmp.bam")
    slice_cmd = "samtools view -b {} {} > {}"
    slice_cmd = slice_cmd.format(bam_file, chroms, tmp)
    system(slice_cmd)
    system("samtools index {}".format(tmp))
    cmd = " bam2hints --intronsonly --in {} --out {}"
    cmd = cmd.format(tmp, f)
    system(cmd)


def cat_bam_2_hints(target, intron_hints_tree, intron_gff_path):
    combined_bam_hints = get_tmp(target, global_dir=True, name="combined_bam_hints")
    combined = catFiles(intron_hints_tree.listFiles(), combined_bam_hints)
    cmd = "cat {} | sort -n -k4,4 | sort -s -n -k5,5 | sort -s -n -k3,3 | sort -s -k1,1 | join_mult_hints.pl > {}"
    cmd = cmd.format(combined_bam_hints, intron_gff_path)
    system(cmd)
    assert os.path.getsize(intron_gff_path) > 10000, cmd


def cat_hints(target, exon_gff_path, intron_gff_path, db_path, genome, genome_fasta):
    hints_file = get_tmp(target, global_dir=True, name="hints_file")
    system("cat {} {} > {}".format(exon_gff_path, intron_gff_path, hints_file))
    target.setFollowOnTargetFn(load_db, args=[hints_file, db_path, genome, genome_fasta])


def load_db(target, hints_file, db_path, genome, genome_fasta):
    cmd = "load2sqlitedb --noIdx --species={} --dbaccess={} {}"
    fa_cmd = cmd.format(genome, db_path, genome_fasta)
    system(fa_cmd)
    hints_cmd = cmd.format(genome, db_path, hints_file)
    system(hints_cmd)


def filter_wrapper(target, bam_paths, db_path, genome, genome_fasta):
    file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), "filter_file_tree"))
    for p in bam_paths:
        target.addChildTargetFn(filter_bam, memory=8 * 1024 ** 3, cpu=2,
                                args=[p, file_tree])
    wig_path = get_tmp(target, global_dir=True, name="wig_path")
    target.setFollowOnTargetFn(merge_bam_to_wig, memory=8 * 1024 ** 3, cpu=2,
                               args=[file_tree, wig_path, db_path, genome, genome_fasta])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--fasta", required=True)
    bamfiles = parser.add_mutually_exclusive_group(required=True)
    bamfiles.add_argument("--bamFiles", nargs="+", help="bamfiles being used", dest="bams")
    bamfiles.add_argument("--bamFofn", help="File containing list of bamfiles", dest="bams")
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    if not isinstance(args.bams, list):
        if not os.path.exists(args.bams):
            raise RuntimeError("ERROR: bamFofn does not exist.")
        args.bams = [x.rstrip() for x in open(args.bams)]
    s = Stack(Target.makeTargetFn(filter_wrapper, memory=8 * 1024 ** 3,
                                  args=[args.bams, args.database, args.genome, args.fasta]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from augustus.build_hints_db import *
    main()
