"""
jobTree wrapper for phastCons. Takes as input a hal file and a reference genome as well as a model file produced by
halPhyloPTrain.py.
Requires you have the hal tools in your path as well as the Phast package.
Will run phastCons twice, following the best practices guide here:
http://compgen.cshl.edu/phast/phastCons-HOWTO.html
Initially, the alignment will be split up into 1mb parts. Then, the conserved/nonconserved models will be estimated
using the first run of phastCons on each alignment. These are then merged using phyloBoot.
"""

import os
import argparse
import itertools
import shutil
import random
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
    parser.add_argument('output_bed', help='Output BED path for most conserved elements.')
    parser.add_argument('output_wig', help='Output wiggle tracks.')
    parser.add_argument('--target-coverage', default='0.05',
                        help='target coverage parameter for phastCons. (default: %(default)s)')
    parser.add_argument('--expected-length', default='20',
                        help='expected length parameter for phastCons. (default: %(default)s)')
    parser.add_argument('--windows', default='1000000,0',
                        help='windows parameter to pass on to msa_view (default: %(default)s)')
    parser.add_argument('--between-blocks', default='5000',
                        help='between blocks parameter to pass on to msa_view (default: %(default)s)')
    parser.add_argument('--min-informative', default='1000',
                        help='min informative parameter to pass on to msa_view (default: %(default)s)')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


def phastcons_pipeline_wrapper(target, args):
    """
    Main pipeline wrapper. Calls out to phast_subset commands, which finally will return to the function passed.
    """
    args.ref_fasta_path = get_ref_genome_fasta(args.hal, args.ref_genome, target.getGlobalTempDir())
    target.setFollowOnTargetFn(extract_maf_wrapper, args=(args, phastcons_estimate_models_wrapper))


def phastcons_estimate_models_wrapper(target, args, split_ss_dict, num_samples=250):
    """
    First pass of phastCons we estimate conserved and nonconserved models. This is done on a random subset of alignments
    in order to save computational effort.
    """
    gc_content = extract_gc_content(args.model)
    output_model_tmp_dir = os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString())
    os.mkdir(output_model_tmp_dir)
    split_ss_list = []
    for chromosome, split_ss_dir in split_ss_dict.iteritems():
        for split_ss in os.listdir(split_ss_dir):
            split_ss_list.append(os.path.join(split_ss_dir, split_ss))
    sampled_ss = random.sample(split_ss_list, k=num_samples)
    for split_ss_path in sampled_ss:
        target.addChildTargetFn(phastcons_estimate_models, args=(args.model, split_ss_path, args.phastcons_options, 
                                                                 gc_content, output_model_tmp_dir))
    target.setFollowOnTargetFn(merge_conserved_nonconserved_models, args=(args, split_ss_dict, model_dir))


def phastcons_estimate_models(target, init_model, split_ss_path, phastcons_options, gc_content, output_model_tmp_dir):
    """
    Run phastcons on randomly selected sub-alignments to produce a conserved and nonconserved model.
    """
    base_cmd = 'phastCons {} --gc {} --estimate-trees {} --no-post-probs {} {}'
    cmd = base_cmd.format(phastcons_options, gc_content, output_model_tmp_dir, split_ss_path, init_model)
    system(cmd)


def merge_conserved_nonconserved_models(target, args, split_ss_dict, model_dir):
    """
    Merged conserved and nonconserved models using phyloBoot
    """
    models = {'cons': os.path.join(target.getGlobalTempDir(), 'ave.cons.mod'),
              'noncons': os.path.join(target.getGlobalTempDir(), 'ave.noncons.mod')}
    for mode, outf in models.iteritems():
        files = [f for f in os.listdir(model_dir) if mode in f]
        fofn = os.path.join(target.getGlobalTempDir(), mode + '.fofn')
        with open(fofn, 'w') as fofn_h:
            for x in files:
                fofn_h.write(x + '\n')
        cmd = 'phyloBoot --read-mods {} --output-average {}'.format(fofn, outf)
        system(cmd)
    models_str = " ".join(models.values())
    target.setFollowOnTargetFn(phastcons_wrapper, args=(args, split_ss_dict, models_str))


def phastcons_wrapper(target, args, split_ss_dict, models_str):
    """
    Wrapper for phastcons function.
    """
    output_bed_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    output_wig_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString()))
    for chromosome, split_ss_dir in split_ss_dict.iteritems():
        for split_ss in os.listdir(split_ss_dir):
            bed_path = output_bed_tree.getTempFile(suffix=split_ss + '.bed')
            wig_path = output_wig_tree.getTempFile(suffix=split_ss + '.wig')
            split_ss_path = os.path.join(split_ss_dir, split_ss)
            target.addChildTargetFn(phastcons, args=(split_ss_path, bed_path, wig_path, models_str, args.phastcons_options))
    target.setFollowOnTargetFn(cat_phastcons, args=(args, output_bed_tree, output_wig_tree))


def phastcons(target, split_ss_path, bed_path, wig_path, models_str, phastcons_options):
    """
    Runs phastcons for the second time, actually producing output.
    """
    cmd = 'phastCons {} --most-conserved {} --score {} {} > {}'.format(phastcons_options, bed_path, split_ss_path, models_str, wig_path)
    system(cmd)


def cat_phastcons(target, args, output_bed_tree, output_wig_tree):
    """
    Concatenates final phastcons output into one gff file.
    """
    for p, t in zip(*[[args.output_bed, args.output_wig]], [output_bed_tree, output_wig_tree]]):
        all_files = t.listFiles()
        fofn = os.path.join(target.getGlobalTempDir(), '{}_fofn'.format(os.path.basename(p)))
        with open(fofn, 'w') as outf:
            for x in all_files:
                outf.write(x + '\n')
        tmp_p = os.path.join(target.getGlobalTempDir(), os.path.basename(p) + '.tmp')
        cat_cmd = 'cat {} | xargs -n 50 cat | sort -k1,1 -k2,2n > {}'.format(fofn, tmp_p)
        system(cat_cmd)
        os.rename(tmp_p, p)


def main():
    args = parse_args()
    args.target_genomes = extract_model_tree(args.model) - set([args.ref_genome])
    args.msa_split_options = " ".join(['--windows', args.windows, '--between-blocks', args.between_blocks, 
                                       '--min-informative', args.min_informative])
    args.phastcons_options = " ".join(['--target-coverage', args.target_coverage, '--expected-length',
                                       args.expected_length])
    s = Stack(Target.makeTargetFn(phastcons_pipeline_wrapper, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from phast.phastcons import *
    main()