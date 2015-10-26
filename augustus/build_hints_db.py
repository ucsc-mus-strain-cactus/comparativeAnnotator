"""
This program takes a series of BAM files and converts them to a Augustus hints database.
"""

import sys
import os
import pysam
import argparse
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


def filter_bam(target, path, file_tree):
    paired = "--paired --pairwiseAlignments" if bam_is_paired(path) is True else ""
    cmd = "samtools sort -O bam -T {} -n {} | filterBam --uniq {} --in /dev/stdin --out {}"
    tmp = os.path.join(target.getLocalTempDir(), getRandomAlphaNumericString(10))
    cmd = cmd.format(tmp, path, paired, file_tree.getTempFile())
    system(cmd)


def merge_bam_to_wig(target, file_tree, wig_path, db_path, genome, genome_fasta):
    bam_files = file_tree.listFiles()
    cmd = "bamtools merge {} | bam2wig /dev/stdin > {}".format(" ".join(bam_files), wig_path)
    system(cmd)
    exon_gff_path = os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString(10))
    intron_gff_path = os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString(10))
    target.addChildTargetFn(wig_2_hints, memory=8 * 1024 ** 3, cpu=2,
                            args=[wig_path, exon_gff_path])
    target.addChildTargetFn(bam_2_hints_wrapper, args=[bam_files, intron_gff_path])
    target.setFollowOnTargetFn(cat_hints, args=[exon_gff_path, intron_gff_path, db_path, genome, genome_fasta])


def wig_2_hints(target, wig_path, exon_gff_path):
    
    cmd = ('cat {} | wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep '
           '--UCSC=/dev/null --radius=4.5 --pri=4 --strand="." > {}')
    cmd = cmd.format(wig_path, exon_gff_path)
    system(cmd)


def bam_2_hints_wrapper(target, bam_files, intron_gff_path):
    intron_hints_tree = TempFileTree(target.getGlobalTempDir())
    for f in bam_files:
        target.addChildTargetFn(bam_2_hints, memory=8 * 1024 ** 3, cpu=1, 
                                args=[f, intron_hints_tree])
    
    target.setFollowOnTargetFn(cat_bam_2_hints, memory=8 * 1024 ** 3, cpu=4,
                                args=[intron_hints_tree, intron_gff_path])


def bam_2_hints(target, bam_file, intron_hints_tree):
    cmd = "bam2hints --intronsonly --in {} --out {}"
    cmd = cmd.format(bam_file, intron_hints_tree.getTempFile())
    system(cmd)


def cat_bam_2_hints(target, intron_hints_tree, intron_gff_path):
    cmd = "cat {} | sort -n -k4,4 | sort -s -n -k5,5 | sort -s -n k3,3 | sort -s -k1,1 | join_mult_hints.pl > {}"
    cmd = cmd.format(" ".join(intron_hints_tree.listFiles()), intron_gff_path)


def cat_hints(target, exon_gff_path, intron_gff_path, db_path, genome, genome_fasta):
    hints_file = os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString(10))
    system("cat {} {} > {}".format(exon_gff_path, intron_gff_path, hints_file))
    target.setFollowOnTargetFn(load_db, args=[hints_file, db_path, genome, genome_fasta])


def load_db(taret, hints_file, db_path, genome, genome_fasta):
    cmd = "load2sqlitedb --noIdx --species={} --dbaccess={} {}"
    fa_cmd = cmd.format(genome, db_path, genome_fasta)
    system(fa_cmd)
    hints_cmd = cmd.format(genome, db_path, hints_file)
    system(hints_cmd)


def filter_wrapper(target, bam_paths, db_path, genome, genome_fasta):
    file_tree = TempFileTree(target.getGlobalTempDir())
    for p in bam_paths:
        target.addChildTargetFn(filter_bam, memory=8 * 1024 ** 3, cpu=2,
                                args=[p, file_tree])
    wig_path = os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString(10))
    target.setFollowOnTargetFn(merge_bam_to_wig, memory=8 * 1024 ** 3, cpu=2,
                               args=[file_tree, wig_path, db_path, genome, genome_fasta])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--bamFiles", nargs="+", required=True)
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    s = Stack(Target.makeTargetFn(filter_wrapper, args=[args.bamFiles, args.database, args.genome, args.fasta]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from augustus.build_hints_db import *
    main()
