"""
Specific to mouse strain project.

Modifies a basic/comprehensive GTF to have the following properties:

1) gene_id -> parent_gene_id
2) transcript_id -> parent_transcript_id
3) add biotype tag
4) creates a strain specific gene/transcript ID in the following format:
MGP_<strain_id>_<G/T/E><8 digit number>.version (which will be 1 for now)
strain IDs are truncated at 10 characters, and have no underscores, but the last character is retained

"""
import argparse
import sys
import os
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from comparativeAnnotator.database_queries import get_transcript_biotype_map


id_template = 'MGP_{strain_id:.9}{strain_suffix}_{tag_type}{unique_id:010d}.1'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help='Input basic/comprehensive GTF', type=argparse.FileType('r'))
    parser.add_argument('genome', help='strain name')
    parser.add_argument('comp_db', help='comparative annotator database')
    parser.add_argument('out_gtf', help='Output GTF', type=argparse.FileType('w'))
    return parser.parse_args()


def main():
    args = parse_args()
    transcript_biotype_map = get_transcript_biotype_map('C57B6J', args.comp_db)
    gene_counter = 1
    transcript_counter = 1
    seen_genes = {}
    seen_transcripts = {}
    genome = args.genome.replace('_', '')[:-1]
    assert genome.isalnum()
    genome_suffix = args.genome[-1]
    for line in args.gtf:
        line = line.rstrip().split('\t')
        attrs = dict([x.split(' ') for x in line[-1].replace('"', '').split('; ')])
        attrs['parent_gene_id'] = attrs['gene_id']
        attrs['parent_transcript_id'] = attrs['transcript_id']
        # we assume protein coding for CGP transcripts
        attrs['biotype'] = transcript_biotype_map.get(attrs['transcript_id'], 'protein_coding')
        if attrs['parent_gene_id'] not in seen_genes:
            seen_genes[attrs['parent_gene_id']] = gene_counter
            gene_counter += 1
        if attrs['parent_transcript_id'] not in seen_transcripts:
            seen_transcripts[attrs['parent_transcript_id']] = transcript_counter
            transcript_counter += 1
        attrs['gene_id'] = id_template.format(strain_id=genome, strain_suffix=genome_suffix, tag_type='G',
                                              unique_id=seen_genes[attrs['parent_gene_id']])
        attrs['transcript_id'] = id_template.format(strain_id=genome, strain_suffix=genome_suffix, tag_type='T',
                                                    unique_id=seen_transcripts[attrs['parent_transcript_id']])
        if 'exon_id' in attrs:
            attrs['exon_id'] = id_template.format(strain_id=genome, strain_suffix=genome_suffix, tag_type='E',
                                                  unique_id=seen_transcripts[attrs['parent_transcript_id']])
            attrs['exon_id'] += '.' + attrs['exon_number']
        # rebuild the line
        attributes = '; '.join([' '.join([x, '"' + y + '"']) for x, y in attrs.iteritems()])
        line[-1] = attributes
        args.out_gtf.write('\t'.join(line) + '\n')


if __name__ == '__main__':
    main()
