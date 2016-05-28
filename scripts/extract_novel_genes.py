import argparse
import sys
sys.path.insert(0, '/hive/users/ifiddes/ihategit/pipeline/submodules/pycbio')
from pycbio.bio.transcripts import get_transcript_dict


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('cgp_consensus')
    return parser.parse_args()


def main():
    args = parse_args()
    tx_dict = get_transcript_dict(args.cgp_consensus)
    for tx in tx_dict.itervalues():
        if 'jg' in tx.name2 or tx.name2.startswith('g'):
            print '\t'.join(map(str, tx.get_gene_pred()))


if __name__ == '__main__':
    main()
