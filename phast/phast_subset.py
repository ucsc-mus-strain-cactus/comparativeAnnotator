"""
Set of helper functions/jobTree targets for taking a HAL file and producing a large number of small SS files
used to predict conservation by dless or phastCons
"""

import os
import argparse
import itertools
import shutil
from ete3 import Tree
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, popenCatch, TempFileTree, getRandomAlphaNumericString


#############
#   Functions
#############

def get_chromosomes(hal, ref_genome):
    """
    Returns a set of chromosomes present in the reference genome.
    """
    sizes = popenCatch('halStats {} --chromSizes {}'.format(hal, ref_genome))
    sizes = sizes.split("\n")[:-1]  # last line is empty newline
    return {x.split()[0] for x in sizes} 


def get_ref_genome_fasta(hal, ref_genome, tmp_dir):
    """
    Extracts the reference fasta from the hal.
    """
    fasta_path = os.path.join(tmp_dir, '{}.fasta'.format(ref_genome))
    system('hal2fasta {} {} > {}'.format(hal, ref_genome, fasta_path))
    return fasta_path


def extract_model_tree(model_file):
    """
    Extracts the tree from the model file.
    """
    # the last line of a model file is the tree
    lines = open(model_file).readlines()
    l = lines[-1].split("TREE: ")[1]
    model_tree = Tree(l, format=1)
    return set(model_tree.get_leaf_names())


def extract_gc_content(model_file):
    """
    Extracts GC content from a model file.
    """
    lines = open(model_file).readlines()
    l = lines[3].split()
    return int(l[2]) + int(l[3])


#################################
#   Alignment Subsetting Pipeline
#################################


def extract_maf(target, hal, chromosome, ref_genome, target_genomes, maf_path):
    """
    Extracts one chromosome from a MAF.
    """
    base_cmd = 'hal2maf --refSequence {} --refGenome {} --targetGenomes {} --noAncestors {} {}'
    cmd = base_cmd.format(chromosome, ref_genome, ','.join(target_genomes), hal, maf_path)
    system(cmd)


def extract_maf_wrapper(target, args, exit_fn):
    """
    Entry point to pipeline.
    Wrapper for extract_maf, runs extract_maf on each chromosome found by get_chromosomes
    exit_fn should be the function called at the end of this chain of commands (dless, phastCons)
    """
    chromosomes = get_chromosomes(args.hal, args.ref_genome)
    maf_file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    maf_dict = {}
    for chromosome in chromosomes:
        maf_path = maf_file_tree.getTempFile(suffix=chromosome)
        maf_dict[chromosome] = maf_path
        target.addChildTargetFn(extract_maf, args=(args.hal, chromosome, args.ref_genome, args.target_genomes, maf_path))
    target.setFollowOnTargetFn(convert_maf_to_ss_wrapper, args=(args, maf_dict, exit_fn))


def convert_maf_to_ss(target, chromosome, maf_path, ss_path):
    """
    Converts a maf to a ss file, using Joel's hacks.
    """
    tmp_maf_path = os.path.join(target.getLocalTempDir(), chromosome + '.noheader.maf')
    system('sed -e 2d {} > {}'.format(maf_path, tmp_maf_path))  # strip out hal comments from maf
    system('msa_view -o SS {} > {}'.format(tmp_maf_path, ss_path))


def convert_maf_to_ss_wrapper(target, args, maf_dict, exit_fn):
    """
    Wrapper for convert_maf_to_ss, runs convert_maf_to_ss on each maf produced by extract_maf 
    """
    ss_file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    ss_dict = {}
    for chromosome, maf_path in maf_dict.iteritems():
        ss_path = ss_file_tree.getTempFile(suffix=chromosome)
        ss_dict[chromosome] = ss_path
        target.addChildTargetFn(convert_maf_to_ss, args=(chromosome, maf_path, ss_path))
    target.setFollowOnTargetFn(split_ss_wrapper, args=(args, ss_dict, exit_fn))


def split_ss(target, chromosome, ss_path, split_ss_dir, ref_fasta_path, msa_split_options):
    """
    Splits a SS file using msa_split.
    Note: this step will verify that the input SS is valid, I.E. it has more than one sequence. One sequence can occur
    when there is nothing aligned. Seems to be common with unassigned contigs.
    """
    nseqs = int(open(ss_path).next().split()[-1])
    if nseqs > 1:
        base_cmd = 'msa_split {} --refseq {} --in-format SS --out-format SS --out-root {} {}'
        root_path = os.path.join(split_ss_dir, chromosome)
        cmd = base_cmd.format(ss_path, ref_fasta_path, root_path, msa_split_options)
        system(cmd)


def split_ss_wrapper(target, args, ss_dict, exit_fn):
    """
    Wrapper for split_ss, which splits a sufficient statistics files by windows and blocks to make manageable parts.
    """
    split_ss_file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    split_ss_dict = {}
    for chromosome, ss_path in ss_dict.iteritems():
        split_ss_dir = split_ss_file_tree.getTempDirectory()
        split_ss_dict[chromosome] = split_ss_dir
        target.addChildTargetFn(split_ss, args=(chromosome, ss_path, split_ss_dir, args.ref_fasta_path, 
                                                args.msa_split_options))
    target.setFollowOnTargetFn(exit_fn, args=(args, split_ss_dict))
