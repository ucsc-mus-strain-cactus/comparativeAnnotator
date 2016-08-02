"""
Specific to mouse strain project.

Modifies a basic/comprehensive GTF to have the following properties:

1) gene_id -> parent_gene_id
2) transcript_id -> parent_transcript_id
3) add biotype tag
4) creates a strain specific gene/transcript ID in the following format:
MGP_<strain_id>_<G/T/E><7 digit number>
strain IDs are truncated at 10 characters, and have no underscores, but the last character is retained

"""
import argparse
import sys
import os
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from comparativeAnnotator.database_queries import get_transcript_biotype_map, get_gene_biotype_map


id_template = 'MGP_{strain_id:.10}_{tag_type}{unique_id:07d}'


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('comp_gtf', help='Input comprehensive (combined) GTF', type=argparse.FileType('r'))
    parser.add_argument('basic_gp', help='Input basic GP', type=argparse.FileType('r'))
    parser.add_argument('comp_gp', help='Input comprehensive (combined) GP', type=argparse.FileType('r'))
    parser.add_argument('genome', help='strain name')
    parser.add_argument('comp_db', help='comparative annotator database')
    parser.add_argument('out_comp_gtf', help='Output comprehensive GTF', type=argparse.FileType('w'))
    parser.add_argument('out_basic_gtf', help='Output basic GTF', type=argparse.FileType('w'))
    parser.add_argument('out_comp_table', help='Output comprehensive data table', type=argparse.FileType('w'))
    parser.add_argument('out_basic_table', help='Output basic data table', type=argparse.FileType('w'))
    parser.add_argument('out_comp_gp', help='Output comprehensive genePred', type=argparse.FileType('w'))
    parser.add_argument('out_basic_gp', help='Output basic genePred', type=argparse.FileType('w'))
    return parser.parse_args()


def evaluate_gtf_line(line, transcript_biotype_map, gene_biotype_map, seen_genes, gene_counter, seen_transcripts,
                      transcript_counter, genome, tx_table):
    """evaluates one GTF line, parsing and adding results to the database"""
    line = line.rstrip().split('\t')
    attrs = dict([x.split(' ') for x in line[-1].replace('"', '').split('; ')])
    attrs['parent_gene_id'] = attrs['gene_id']
    attrs['parent_transcript_id'] = attrs['transcript_id']
    # we assume protein coding for CGP transcripts
    attrs['transcript_biotype'] = transcript_biotype_map.get(attrs['transcript_id'], 'unknown_likely_coding')
    attrs['gene_biotype'] = gene_biotype_map.get(attrs['gene_id'], 'unknown_likely_coding')
    if attrs['parent_gene_id'] not in seen_genes:
        seen_genes[attrs['parent_gene_id']] = gene_counter
        gene_counter += 1
    if attrs['parent_transcript_id'] not in seen_transcripts:
        seen_transcripts[attrs['parent_transcript_id']] = transcript_counter
        transcript_counter += 1
    attrs['gene_id'] = id_template.format(strain_id=genome, tag_type='G',
                                          unique_id=seen_genes[attrs['parent_gene_id']])
    attrs['transcript_id'] = id_template.format(strain_id=genome, tag_type='T',
                                                unique_id=seen_transcripts[attrs['parent_transcript_id']])
    if 'exon_id' in attrs:
        attrs['exon_id'] = id_template.format(strain_id=genome, tag_type='E',
                                              unique_id=seen_transcripts[attrs['parent_transcript_id']])
        attrs['exon_id'] += '.' + attrs['exon_number']
    # rebuild the line
    attributes = '; '.join([' '.join([x, '"' + y + '"']) for x, y in attrs.iteritems()]) + ';'
    line[-1] = attributes
    if line[2] == 'transcript':
        tx_table[attrs['transcript_id']] = [attrs['transcript_id'], attrs['parent_transcript_id'],
                                            attrs['gene_id'], attrs['parent_gene_id'],
                                            attrs['gene_biotype'], attrs['transcript_biotype']]
    return gene_counter, transcript_counter, line, attrs


def rename_genepreds(tx_table, gp_recs, out_gp_handle):
    """rename the genePred using the tx_table"""
    # generate a mapping of original names to new names
    original_name_map = {parent_transcript_id: (transcript_id, gene_id) for transcript_id, parent_transcript_id,
                         gene_id, parent_gene_id, gene_biotype, transcript_biotype in tx_table.itervalues()}
    for rec in gp_recs:
        rec = rec.split()
        new_name, new_name2 = original_name_map[rec[0]]
        rec[0] = new_name
        rec[11] = new_name2
        out_gp_handle.write('\t'.join(rec) + '\n')


def main():
    args = parse_args()
    transcript_biotype_map = get_transcript_biotype_map('C57B6J', args.comp_db)
    gene_biotype_map = get_gene_biotype_map('C57B6J', args.comp_db)
    gene_counter = 1
    transcript_counter = 1
    seen_genes = {}
    seen_transcripts = {}
    genome = args.genome.replace('_', '')
    assert genome.isalnum(), genome
    tx_table = {}
    basic_gp_recs = list(args.basic_gp)
    comp_gp_recs = list(args.comp_gp)
    basic_ids = {x.split()[0] for x in basic_gp_recs}
    for line in args.comp_gtf:
        gene_counter, transcript_counter, updated_line, attrs = evaluate_gtf_line(line, transcript_biotype_map,
                                                                                  gene_biotype_map, seen_genes,
                                                                                  gene_counter, seen_transcripts,
                                                                                  transcript_counter, genome, tx_table)
        args.out_comp_gtf.write('\t'.join(updated_line) + '\n')
        if attrs['parent_transcript_id'] in basic_ids:
            args.out_basic_gtf.write('\t'.join(updated_line) + '\n')
    args.out_comp_table.write('\t'.join(['#transcript_id', 'parent_transcript_id', 'gene_id', 'parent_gene_id',
                                         'gene_biotype', 'transcript_biotype']) + '\n')
    args.out_basic_table.write('\t'.join(['#transcript_id', 'parent_transcript_id', 'gene_id', 'parent_gene_id',
                                          'gene_biotype', 'transcript_biotype']) + '\n')
    for rec in tx_table.itervalues():
        args.out_comp_table.write('\t'.join(rec) + '\n')
        if rec[1] in basic_ids:
            args.out_basic_table.write('\t'.join(rec) + '\n')
    rename_genepreds(tx_table, basic_gp_recs, args.out_basic_gp)
    rename_genepreds(tx_table, comp_gp_recs, args.out_comp_gp)


if __name__ == '__main__':
    main()
