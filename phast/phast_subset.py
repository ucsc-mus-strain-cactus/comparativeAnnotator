"""
Set of helper functions/jobTree targets for taking a HAL file and producing a large number of small SS files
used to predict conservation by dless or phastCons
"""

import os
import argparse
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system
from phast.phast_functions import *
from lib.general_lib import mkdir_p


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hal', help='HAL alignment file.')
    parser.add_argument('ref_genome', help='Reference genome.')
    parser.add_argument('model', help='Model file produced by phyloFit/halPhyloPTrain.py.')
    parser.add_argument('output_dir', help='Location to write the split up files to.')
    parser.add_argument('--ref-fasta-path', default=None,
                        help='Path to reference genome FASTA. If not provided, it will be extracted from the HAL.')
    parser.add_argument('--windows', default='1000000,0',
                        help='windows parameter to pass on to msa_view (default: %(default)s)')
    parser.add_argument('--between-blocks', default='5000',
                        help='between blocks parameter to pass on to msa_view (default: %(default)s)')
    parser.add_argument('--min-informative', default='1000',
                        help='min informative parameter to pass on to msa_view (default: %(default)s)')
    Stack.addJobTreeOptions(parser)
    return parser.parse_args()


def extract_maf(target, hal, chromosome, ref_genome, target_genomes, maf_path):
    """
    Extracts one chromosome from a MAF.
    """
    base_cmd = 'hal2maf --refSequence {} --refGenome {} --targetGenomes {} --noAncestors {} {}'
    cmd = base_cmd.format(chromosome, ref_genome, ','.join(target_genomes), hal, maf_path)
    system(cmd)


def subset_hal_pipeline(target, args, out_path=None):
    """
    Entry point to pipeline.
    Wrapper for extract_maf, runs extract_maf on each chromosome found by get_chromosomes
    Other programs that start here can adjust the output path so its different than their own.
    """
    if out_path is not None:
        args.output_dir = out_path
    chromosomes = get_chromosomes(args.hal, args.ref_genome)
    maf_file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), 'extracted_mafs'))
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
    ss_file_tree = TempFileTree(os.path.join(target.getGlobalTempDir(), 'maf_to_ss'))
    ss_dict = {}
    for chromosome, maf_path in maf_dict.iteritems():
        ss_path = ss_file_tree.getTempFile(suffix=chromosome)
        ss_dict[chromosome] = ss_path
        target.addChildTargetFn(convert_maf_to_ss, args=(chromosome, maf_path, ss_path))
    target.setFollowOnTargetFn(split_ss_wrapper, args=(args, ss_dict))


def split_ss(target, chromosome, ss_path, out_dir, ref_fasta_path, msa_split_options):
    """
    Splits a SS file using msa_split.
    Note: this step will verify that the input SS is valid, I.E. it has more than one sequence. One sequence can occur
    when there is nothing aligned. Seems to be common with unassigned contigs.
    """
    nseqs = extract_nseqs(ss_path)
    if nseqs > 1:
        base_cmd = 'msa_split {} --refseq {} --in-format SS --out-format SS --out-root {} {}'
        cmd = base_cmd.format(ss_path, ref_fasta_path, out_dir, msa_split_options)
        system(cmd)


def split_ss_wrapper(target, args, ss_dict):
    """
    Wrapper for split_ss, which splits a sufficient statistics files by windows and blocks to make manageable parts.
    """
    if args.ref_fasta_path is None:
        args.ref_fasta_path = get_ref_genome_fasta(args.hal, args.ref_genome, target.getGlobalTempDir())
    for chromosome, ss_path in ss_dict.iteritems():
        out_dir = os.path.join(args.output_dir, chromosome) + '/'  # need to add trailing slash for msa_split
        mkdir_p(out_dir)
        target.addChildTargetFn(split_ss, args=(chromosome, ss_path, out_dir, args.ref_fasta_path,
                                                args.msa_split_options))


def main():
    args = parse_args()
    args.target_genomes = extract_model_tree(args.model) - set([args.ref_genome])
    args.msa_split_options = ' '.join(['--windows', args.windows, '--between-blocks', args.between_blocks,
                                       '--min-informative', args.min_informative])
    s = Stack(Target.makeTargetFn(subset_hal_pipeline, args=(args,)))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from phast.phast_subset import *
    main()
