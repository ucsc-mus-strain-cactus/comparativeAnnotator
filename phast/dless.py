"""
jobTree wrapper for dless. Takes as input a hal file and a reference genome as well as a model file produced by
halPhyloPTrain.py.
Requires you have the hal tools in your path as well as the Phast package.
"""

import os
import argparse
import itertools
import shutil
from ete3 import Tree
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, popenCatch, TempFileTree, getRandomAlphaNumericString


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hal', help='HAL alignment file.')
    parser.add_argument('ref_genome', help='Reference genome.')
    parser.add_argument('model', help='Model file produced by phyloFit/halPhyloPTrain.py.')
    parser.add_argument('output_gff', help='Output GFF path.')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


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


def dless_pipeline_wrapper(target, args):
    args.ref_fasta_path = get_ref_genome_fasta(args.hal, args.ref_genome, target.getGlobalTempDir())
    target.setFollowOnTargetFn(extract_maf_wrapper, args=(args,))


def extract_maf(target, hal, chromosome, ref_genome, target_genomes, maf_path):
    """
    Extracts one chromosome from a MAF.
    """
    base_cmd = 'hal2maf --refSequence {} --refGenome {} --targetGenomes {} --noAncestors {} {}'
    cmd = base_cmd.format(chromosome, ref_genome, ','.join(target_genomes), hal, maf_path)
    system(cmd)


def extract_maf_wrapper(target, args):
    """
    Wrapper for extract_maf, runs extract_maf on each chromosome found by get_chromosomes
    """
    chromosomes = get_chromosomes(args.hal, args.ref_genome)
    maf_file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    maf_dict = {}
    for chromosome in chromosomes:
        maf_path = maf_file_tree.getTempFile(suffix=chromosome)
        maf_dict[chromosome] = maf_path
        target.addChildTargetFn(extract_maf, args=(args.hal, chromosome, args.ref_genome, args.target_genomes, maf_path))
    target.setFollowOnTargetFn(convert_maf_to_ss_wrapper, args=(args, maf_dict))


def convert_maf_to_ss(target, chromosome, maf_path, ss_path):
    """
    Converts a maf to a ss file, using Joel's hacks.
    """
    tmp_maf_path = os.path.join(target.getLocalTempDir(), chromosome + '.noheader.maf')
    system('sed -e 2d {} > {}'.format(maf_path, tmp_maf_path))  # strip out hal comments from maf
    system('msa_view -o SS {} > {}'.format(tmp_maf_path, ss_path))


def convert_maf_to_ss_wrapper(target, args, maf_dict):
    """
    Wrapper for convert_maf_to_ss, runs convert_maf_to_ss on each maf produced by extract_maf 
    """
    ss_file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    ss_dict = {}
    for chromosome, maf_path in maf_dict.iteritems():
        ss_path = ss_file_tree.getTempFile(suffix=chromosome)
        ss_dict[chromosome] = ss_path
        target.addChildTargetFn(convert_maf_to_ss, args=(chromosome, maf_path, ss_path))
    target.setFollowOnTargetFn(split_ss_wrapper, args=(args, ss_dict))


def split_ss(target, chromosome, ss_path, split_ss_dir, ref_fasta_path):
    """
    Splits a SS file using msa_split.
    """
    base_cmd = ('msa_split {} --refseq {} --in-format SS --out-format SS --windows 1000000,1000 --between-blocks 5000'
                ' --out-root {}')
    root_path = os.path.join(split_ss_dir, chromosome)
    cmd = base_cmd.format(ss_path, ref_fasta_path, root_path)
    system(cmd)


def split_ss_wrapper(target, args, ss_dict):
    """
    Wrapper for split_ss, which splits a sufficient statistics files by windows and blocks to make manageable parts.
    """
    split_ss_file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    split_ss_dict = {}
    for chromosome, ss_path in ss_dict.iteritems():
        split_ss_dir = split_ss_file_tree.getTempDirectory()
        split_ss_dict[chromosome] = split_ss_dir
        target.addChildTargetFn(split_ss, args=(chromosome, ss_path, split_ss_dir, args.ref_fasta_path))
    target.setFollowOnTargetFn(dless_wrapper, args=(args, split_ss_dict))


def dless(target, split_ss_path, gff_path, model):
    """
    Main function for running dless. Strips all headers out of final gff.
    """
    system('dless {} {} | sed "/^#/ d" {}'.format(split_ss_path, model, gff_path))


def dless_wrapper(target, args, split_ss_dict):
    """
    Wrapper for dless function.
    """
    output_gff_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    for chromosome, split_ss_dir in split_ss_dict.iteritems():
        for split_ss in os.listdir(split_ss_dir):
            gff_path = output_gff_tree.getTempFile(suffix=split_ss + '.gff')
            split_ss_path = os.path.join(split_ss_dir, split_ss)
            target.addChildTargetFn(dless, args=(split_ss_path, gff_path, args.model))
    target.setFollowOnTargetFn(cat_dless, args=(args, output_gff_tree))


def cat_dless(target, args, output_gff_tree):
    """
    Concatenates final dless output into one gff file.
    """
    all_gffs = output_gff_tree.listFiles()
    gff_fofn = os.path.join(target.getGlobalTempDir(), 'gff_fofn')
    with open(gff_fofn, 'w') as outf:
        for x in all_gffs:
            outf.write(x + '\n')
    concat_tmp = os.path.join(target.getGlobalTempDir(), 'concat_regions.gff')
    cat_cmd = 'cat {} | xargs -n 50 >> {}'.format(gff_fofn, concat_tmp)
    system(cat_cmd)
    assert False


def main():
    args = parse_args()
    args.target_genomes = extract_model_tree(args.model) - set([args.ref_genome])
    s = Stack(Target.makeTargetFn(dless_pipeline_wrapper, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from phast.dless import *
    main()