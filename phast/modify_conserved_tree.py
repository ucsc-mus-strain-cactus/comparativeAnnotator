"""
Simple script that modifies a phast model file to adjust all branch lengths in a set of genomes to 1
"""
import sys
import os
import argparse
from ete3 import Tree

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomes', nargs='+', help='Genomes to change branch lengths. Should not be paraphyletic',
                        default=['C57BL_6NJ', 'NZO_HlLtJ', 'NOD_ShiLtJ', 'FVB_NJ', 'LP_J', '129S1_SvImJ', 'AKR_J', 'BALB_cJ', 'A_J', 'DBA_2J', 'CBA_J', 'C3H_HeJ', 'WSB_EiJ', 'CAST_EiJ', 'PWK_PhJ', 'SPRET_EiJ', 'CAROLI_EiJ', 'Pahari_EiJ'])
    parser.add_argument('model', type=argparse.FileType('r'))
    return parser.parse_args()


def main():
    args = parse_args()
    lines = args.model.readlines()
    t = Tree(lines[-1].split('TREE: ')[1], format=1)
    # isolate only the genomes we have this time
    genomes = [x for x in args.genomes if x in t.get_leaf_names()]
    assert len(genomes) > 0
    # move to the ancestor of all of these, set all children distances to 1
    anc = t.get_common_ancestor(genomes)
    anc._set_dist(1)
    for node in anc.get_descendants():
        node._set_dist(1)
    # rebuild the new model
    for l in lines[:-1]:
        sys.stdout.write(l)
    new_tree = 'TREE: ' + t.write()
    print new_tree


if __name__ == '__main__':
    main()
