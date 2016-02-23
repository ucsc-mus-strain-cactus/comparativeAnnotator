"""
Useful functions for working with HAL formatted files for phastCons/dless.
"""
import os
from ete3 import Tree
from sonLib.bioio import system, popenCatch, TempFileTree, getRandomAlphaNumericString


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
    return float(l[2]) + float(l[3])


def extract_nseqs(ss_file):
    """
    Extracts nseqs from sufficient statistics file.
    """
    nseqs = int(open(ss_file).next().rstrip().split(' = ')[-1])
    return nseqs


def read_subalignment_dir(path):
    """
    Builds a dict from a SS subalignment dir produced by phast_subset.
    Structure: {ref_chrom: (chrom, start, stop): path}
    """
    r = {}
    for base_dir, dirs, files in os.walk(os.path.abspath(path)):
        if files:
            for split_ss in files:
                chrom, start_stop = split_ss.rsplit(".", 2)[:-1]  # some chrom names have a period
                start, stop = map(int, start_stop.split("-"))
                r[(chrom, start, stop)] = os.path.join(base_dir, split_ss)
    return r
