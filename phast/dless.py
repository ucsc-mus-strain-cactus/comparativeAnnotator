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
from phast.phast_subset import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hal', help='HAL alignment file.')
    parser.add_argument('ref_genome', help='Reference genome.')
    parser.add_argument('model', help='Model file produced by phyloFit/halPhyloPTrain.py.')
    parser.add_argument('output_gff', help='Output GFF path.')
    parser.add_argument('--windows', default='1000000,1000',
                        help='windows parameter to pass on to msa_view (default: %(default)s)')
    parser.add_argument('--between-blocks', default='5000',
                        help='between blocks parameter to pass on to msa_view (default: %(default)s)')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


def dless_pipeline_wrapper(target, args):
    """
    Main pipeline wrapper. Calls out to phast_subset commands, which finally will return to the function passed.
    """
    args.ref_fasta_path = get_ref_genome_fasta(args.hal, args.ref_genome, target.getGlobalTempDir())
    target.setFollowOnTargetFn(extract_maf_wrapper, args=(args, dless_wrapper))


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


def dless(target, split_ss_path, gff_path, model):
    """
    Main function for running dless. Strips all headers out of final gff.
    """
    system('dless {} {} | sed "/^#/ d" > {}'.format(split_ss_path, model, gff_path))


def cat_dless(target, args, output_gff_tree):
    """
    Concatenates final dless output into one gff file.
    """
    all_gffs = output_gff_tree.listFiles()
    gff_fofn = os.path.join(target.getGlobalTempDir(), 'gff_fofn')
    with open(gff_fofn, 'w') as outf:
        for x in all_gffs:
            outf.write(x + '\n')
    cat_cmd = 'cat {} | xargs -n 50 cat >> {}'.format(gff_fofn, args.output_gff)
    system(cat_cmd)


def main():
    args = parse_args()
    args.target_genomes = extract_model_tree(args.model) - set([args.ref_genome])
    args.msa_split_options = " ".join(['--windows', args.windows, '--between-blocks', args.between_blocks])
    s = Stack(Target.makeTargetFn(dless_pipeline_wrapper, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from phast.dless import *
    main()