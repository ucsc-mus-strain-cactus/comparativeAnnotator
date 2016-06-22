"""
Generates a basic/comprehensive split for a gene set. This is adapted from the following perl script:

https://github.com/Ensembl/ensembl-analysis/blob/master/scripts/Merge/label_gencode_basic_transcripts.pl

The notes from that script are below, and lay out the decision tree process used here.

In order to adapt this, I need to know if a given transcript is incomplete. The genePred format does not tell me
if the tss/tts are complete or not. To find this information, I have different sources:
AugustusTMR/CGP: Is present in the GTF by the presence/absence of a tss/tts feature
transMap: Is present in the source GTF, and can use the PSL to determine if the mapped over version is the same.


1) Protein-coding transcripts

a)We have preference to these transcript biotypes :

"IG_D_gene" ,"IG_J_gene" ,"IG_C_gene", "IG_V_gene" ,"protein_coding",
"TR_C_gene","TR_J_gene","TR_V_gene" ,"TR_D_gene", "polymorphic_pseudogene"

 If we have these biotypes, in full-length transcripts we put them ALL
in the basic set.


b)If we don't get any of these,    we check these biotypes:

 "IG_D_gene" , "IG_J_gene" , "IG_C_gene" , "IG_V_gene"
,"protein_coding" , "TR_C_gene","TR_J_gene", "TR_V_gene" ,"TR_D_gene" ,
"non_stop_decay" ,"nonsense_mediated_decay" , "polymorphic_pseudogene"

to get ONE.  We want to put the longest CDS in the basic set.
if there is more than one with longest CDS we put the biggest length in the basic set.
if there is more than one with biggest length we put the longest genomic
length in the basic set.(if there is more than 1 of these too, I get the first I find.)



2) Non - coding transcripts

a)We have preference to these transcript biotypes :
 "antisense", "Mt_tRNA","Mt_rRNA" , "rRNA" ,"snoRNA" , "snRNA", "miRNA"
 If we have these biotypes, in full-length transcripts we put them ALL
in the basic set.


b)If we don't get any of these,  we check these biotypes:

"antisense", "Mt_tRNA", "Mt_rRNA","rRNA" , "snoRNA" ,"snRNA",
"processed_transcript", "lincRNA" ,
"3prime_overlapping_ncrna","non_coding" ,
"sense_intronic" ,"sense_overlapping", "miRNA" , "misc_RNA"
to get ONE.  We want to put the one with the biggest length in the basic set.
if there is more than one with biggest length we put the transcript with
the longest genomic length in the basic set.(if there is more than 1 of these too, I get the first I find.)

We want at least 1 transcript from both categories(basic coding and basic non-coding) if there is any
transcript falling in these groups.

3) If we don't get anything from both 1 and 2, last step is that I check the the UCSC the problem category which are the transcripts with these biotypes:

 "retained_intron", "TEC" , "ambiguous_orf" , "disrupted_domain"

I get 1 transcript which will have the biggest length.
If there is more than 1 with the biggest length I get the one with the biggest genomic length! (if there is more than 1 of these too, I get the first I find.)

"""
import argparse
import os
import sys
from collections import defaultdict
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from comparativeAnnotator.database_queries import get_transcript_biotype_map
from pycbio.bio.transcripts import get_transcript_dict
from pycbio.bio.psl import get_alignment_dict
from pycbio.sys.mathOps import format_ratio
from pycbio.sys.fileOps import ensureFileDir
from comparativeAnnotator.comp_lib.gff3 import Gff3Parser
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers

known_biotypes = {
                     'protein_coding'                     : 'coding',
                     'polymorphic_pseudogene'             : 'coding',
                     'IG_D_gene'                          : 'coding',
                     'IG_J_gene'                          : 'coding',
                     'IG_C_gene'                          : 'coding',
                     'IG_V_gene'                          : 'coding',
                     'IG_LV_gene'                         : 'coding',
                     'TR_C_gene'                          : 'coding',
                     'TR_J_gene'                          : 'coding',
                     'TR_V_gene'                          : 'coding',
                     'TR_D_gene'                          : 'coding',
                     'non_stop_decay'                     : 'coding_second_choice',
                     'nonsense_mediated_decay'            : 'coding_second_choice',
                     'rRNA'                               : 'noncoding',
                     'known_ncrna'                        : 'noncoding',
                     'snoRNA'                             : 'noncoding',
                     'snRNA'                              : 'noncoding',
                     'miRNA'                              : 'noncoding',
                     'antisense'                          : 'noncoding',
                     'Mt_tRNA'                            : 'noncoding',
                     'Mt_rRNA'                            : 'noncoding',
                     'sense_intronic'                     : 'noncoding',
                     'sense_overlapping'                  : 'noncoding',
                     'lincRNA'                            : 'noncoding',
                     'macro_lncRNA'                       : 'noncoding_second_choice',
                     'ribozyme'                           : 'noncoding_second_choice',
                     'scaRNA'                             : 'noncoding_second_choice',
                     'scRNA'                              : 'noncoding_second_choice',
                     'sRNA'                               : 'noncoding_second_choice',
                     'vaultRNA'                           : 'noncoding_second_choice',
                     'processed_transcript'               : 'noncoding_second_choice',
                     'misc_RNA'                           : 'noncoding_second_choice',
                     '3prime_overlapping_ncrna'           : 'noncoding_second_choice',
                     'non_coding'                         : 'noncoding_second_choice',
                     'bidirectional_promoter_lncrna'      : 'noncoding_second_choice',
                     'transcribed_processed_pseudogene'   : 'pseudogene_transcribed',
                     'transcribed_unitary_pseudogene'     : 'pseudogene_transcribed',
                     'transcribed_unprocessed_pseudogene' : 'pseudogene_transcribed',
                     'pseudogene'                         : 'pseudogene',
                     'processed_pseudogene'               : 'pseudogene',
                     'unprocessed_pseudogene'             : 'pseudogene',
                     'translated_processed_pseudogene'    : 'pseudogene',
                     'translated_unprocessed_pseudogene'  : 'pseudogene',
                     'unitary_pseudogene'                 : 'pseudogene',
                     'IG_pseudogene'                      : 'pseudogene',
                     'IG_C_pseudogene'                    : 'pseudogene',
                     'IG_D_pseudogene'                    : 'pseudogene',
                     'IG_J_pseudogene'                    : 'pseudogene',
                     'IG_pseudogene'                      : 'pseudogene',
                     'IG_V_pseudogene'                    : 'pseudogene',
                     'TR_J_pseudogene'                    : 'pseudogene',
                     'TR_V_pseudogene'                    : 'pseudogene',
                     'Mt_rRNA_pseudogene'                 : 'pseudogene',
                     'miRNA_pseudogene'                   : 'pseudogene',
                     'misc_RNA_pseudogene'                : 'pseudogene',
                     'rRNA_pseudogene'                    : 'pseudogene',
                     'scRNA_pseudogene'                   : 'pseudogene',
                     'snRNA_pseudogene'                   : 'pseudogene',
                     'snoRNA_pseudogene'                  : 'pseudogene',
                     'tRNA_pseudogene'                    : 'pseudogene',
                     'retained_intron'                    : 'problem',
                     'TEC'                                : 'problem',
                     'ambiguous_orf'                      : 'problem',
                     'disrupted_domain'                   : 'problem',
                     'LRG_gene'                           : 'do_not_use',
                     }


def is_original_cds_stop(tm_tx, ref_tx, tm_psl):
    """
    Does this transMap transcript have its original CDS stop?
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param ref_tx: GenePredTranscript object for reference transcript
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :return: boolean
    """
    for i in xrange(tm_tx.cds_size - 4, tm_tx.cds_size - 1):
        p = tm_tx.chromosome_coordinate_to_cds(tm_psl.query_coordinate_to_target(ref_tx.cds_coordinate_to_transcript(i)))
        if p is None:
            return False
    return True


def is_original_cds_start(tm_tx, ref_tx, tm_psl):
    """
    Does this transMap transcript have its original CDS start?
    :param tm_tx: GenePredTranscript object for transMap transcript
    :param ref_tx: GenePredTranscript object for reference transcript
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :return: boolean
    """
    for i in xrange(3):
        p = tm_tx.chromosome_coordinate_to_cds(tm_psl.query_coordinate_to_target(ref_tx.cds_coordinate_to_transcript(i)))
        if p is None:
            return False
    return True


def is_original_tss(tm_psl, fuzz_distance=25):
    """
    Is this transMap alignment within fuzz_distance from the original transcription start site?
    This can only be calculated inwards from the original start
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :param fuzz_distance: integer distance
    :return: boolean
    """
    if tm_psl.strand == '+':
        return tm_psl.q_start < fuzz_distance
    else:
        return tm_psl.q_size - tm_psl.q_end < fuzz_distance


def is_original_tts(tm_psl, fuzz_distance=25):
    """
    Is this transMap alignment within fuzz_distance from the original transcription termination point?
    This can only be calculated inwards from the original stop
    :param tm_psl: PslRow object for the relationship between tm_tx and ref_tx
    :param fuzz_distance: integer distance
    :return: boolean
    """
    if tm_psl.strand == '-':
        return tm_psl.q_start < fuzz_distance
    else:
        return tm_psl.q_size - tm_psl.q_end < fuzz_distance


def find_complete_transmap(tm_dict, psl_dict, ref_tm_dict, transcript_biotype_map):
    """returns the list of transMap IDs which are complete"""
    r = set()
    for aln_id, tx in tm_dict.iteritems():
        psl = psl_dict[aln_id]
        ref_tx = ref_tm_dict[strip_alignment_numbers(aln_id)]
        tx_biotype = known_biotypes[transcript_biotype_map[ref_tx.name]]
        if 'noncoding' in tx_biotype:
            if is_original_tss(psl) and is_original_tts(psl):
                r.add(aln_id)
        elif 'coding' in tx_biotype:
            if all((is_original_cds_start(tx, ref_tx, psl), is_original_cds_stop(tx, ref_tx, psl),
                   is_original_tss(psl), is_original_tts(psl))):
                r.add(aln_id)
        else:
            r.add(aln_id)
    return r


def extract_reference_completeness(gff3_path, tags=('mRNA_end_NF', 'mRNA_Start_NF', 'cds_end_NF', 'cds_start_NF')):
    """create a set of transcripts which were complete in the reference"""
    ref_gff3 = Gff3Parser(gff3_path)
    ref_gff3 = ref_gff3.parse()
    r = set()
    for gene in ref_gff3.roots:
        for tx in gene.children:
            tx_id = tx.attributes['transcript_id'][0]
            if 'tag' not in tx.attributes:
                r.add(tx_id)
            elif not any(x in tx.attributes['tag'] for x in tags):
                r.add(tx_id)
    return r


def extract_aug_completeness(gtf_file_path, features={'tts', 'tss', 'stop_codon', 'start_codon'}):
    """creates a set of transcripts which are complete in this augustus run. Also returns all IDs found"""
    r = defaultdict(set)
    for l in open(gtf_file_path):
        l = l.split('\t')
        if l[2] in features:
            attributes = l[-1].replace('"', '').replace(';', '').split(' ')
            attr_dict = dict(zip(*[iter(attributes)] * 2))
            r[attr_dict['transcript_id']].add(l[2])
    return {tx_id for tx_id, vals in r.iteritems() if vals == features}, r.viewkeys()


def get_gene_dict(consensus_dict):
    """Create a new dict mapping all transcripts to their parent gene."""
    r = defaultdict(list)
    for tx in consensus_dict.itervalues():
        r[tx.name2].append(tx)
    return r


def get_complete_transcripts(tx_list, complete_ids):
    """returns the subset of these transcripts which are all complete"""
    r = []
    for tx in tx_list:
        if tx.name in complete_ids and strip_alignment_numbers(tx.name) in complete_ids:
            r.append(tx)
    return r


def get_gene_exons(tx_list):
    """get a set of all exons for this gene"""
    r = set()
    for tx in tx_list:
        r.update(tx.exon_intervals)
    return r


def get_coverage_transcripts(tx_list, coverage=80.0):
    """
    Returns the longest transcripts covering the maximum number of exons until coverage is hit
    """
    gene_intervals = get_gene_exons(tx_list)
    tot_len = sum(len(i) for i in gene_intervals)
    sorted_tx_list = sorted(tx_list, key=len, reverse=True)
    r = []
    covered_exons = set()
    # this is a ugly n^3 algorithm
    for tx in sorted_tx_list:
        tot_cov = sum(len(i) for i in covered_exons) if len(covered_exons) > 0 else 0
        if format_ratio(tot_cov, tot_len) >= coverage:
            return r
        r.append(tx)
        for exon in tx.exon_intervals:
            for gene_exon in gene_intervals:
                i = exon.intersection(gene_exon)
                if i is not None:
                    covered_exons.add(i)
    return r


def get_basic_coding(tx_list, complete_ids, transcript_biotype_map):
    """
    Returns all basic coding. If there are any complete coding transcripts, return all. Otherwise, return the one
    longest based on CDS length from the expanded pool of biotypes.
    """
    complete_coding = []
    for tx in get_complete_transcripts(tx_list, complete_ids):
        if known_biotypes[transcript_biotype_map[tx.name]] == 'coding':
            complete_coding.append(tx)
    if len(complete_coding) > 0:
        return complete_coding
    else:
        biotype_list = ['coding', 'coding_second_choice']
        coding_txs = [tx for tx in tx_list if known_biotypes[transcript_biotype_map[tx.name]] in biotype_list]
        if len(coding_txs) == 0:
            return []
        longest = sorted(coding_txs, key=lambda x: x.cds_size, reverse=True)[0]
        return [longest]


def get_basic_noncoding(tx_list, complete_ids, transcript_biotype_map):
    """
    Returns all basic non-coding. If there are any complete non-coding transcripts, return all. Otherwise,
    return the set from the extended pool that together cover at least coverage of the gene exons.
    """
    complete_noncoding = []
    for tx in get_complete_transcripts(tx_list, complete_ids):
        if known_biotypes[transcript_biotype_map[tx.name]] == 'noncoding':
            complete_noncoding.append(tx)
    if len(complete_noncoding) > 0:
        return complete_noncoding
    else:
        tx_subset = [tx for tx in tx_list if known_biotypes[transcript_biotype_map[tx.name]] == 'noncoding']
        if len(tx_subset) > 0:
            return get_coverage_transcripts(tx_subset)
        biotype_list = ['noncoding', 'noncoding_second_choice']
        noncoding_txs = [tx for tx in tx_list if known_biotypes[transcript_biotype_map[tx.name]] in biotype_list]
        if len(noncoding_txs) == 0:
            return []
        longest = sorted(noncoding_txs, key=len, reverse=True)[0]
        return [longest]


def get_basic_pseudogenes(tx_list, transcript_biotype_map):
    """Returns all basic pseudogenes"""
    transcribed_pseudogenes = [tx for tx in tx_list if known_biotypes[transcript_biotype_map[tx.name]] == 'pseudogene_transcribed']
    if len(transcribed_pseudogenes) > 0:
        return transcribed_pseudogenes
    else:
        pseudogenes = [tx for tx in tx_list if 'pseudogene' in known_biotypes[transcript_biotype_map[tx.name]]]
        if len(pseudogenes) == 0:
            return []
        assert len(pseudogenes) == 1
        return pseudogenes[0]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', required=True)
    parser.add_argument('--refGenome', default='C57B6J')
    parser.add_argument('--refGff3', required=True, help='Gencode gff3')
    parser.add_argument('--refGp', required=True, help='Gencode genePred')
    parser.add_argument('--augustusGtf', required=True, help='Augustus TMR GTF')
    #parser.add_argument('--cgpGtf', required=True, help='CGP GTF')
    parser.add_argument('--cgpGp', required=True)
    parser.add_argument('--pacbioGtf', help='pacbio GTF if used')
    parser.add_argument('--transMapGp', required=True, help='TM GP')
    parser.add_argument('--transMapPsl', required=True, help='TM PSL')
    parser.add_argument('--consensus', required=True, help='consensus genePred to divide')
    parser.add_argument('--outBasic', required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    transcript_biotype_map = get_transcript_biotype_map(args.refGenome, args.db)
    tm_dict = get_transcript_dict(args.transMapGp)
    psl_dict = get_alignment_dict(args.transMapPsl)
    ref_tm_dict = get_transcript_dict(args.refGp)
    complete_tm_ids = find_complete_transmap(tm_dict, psl_dict, ref_tm_dict, transcript_biotype_map)
    complete_aug_ids, all_aug_ids = extract_aug_completeness(args.augustusGtf)
    #complete_cgp_ids = extract_aug_completeness(args.cgpGtf)  # we don't have CGP GTF
    complete_cgp_ids = set(get_transcript_dict(args.cgpGp).keys())
    for tx_id in complete_cgp_ids:
        transcript_biotype_map[tx_id] = 'protein_coding'
    if args.pacbioGtf is not None:
        complete_pacbio_ids, all_pacbio_ids = extract_aug_completeness(args.pacbioGtf)
        complete_cgp_ids.update(complete_pacbio_ids)
        for tx_id in all_pacbio_ids:
            transcript_biotype_map[tx_id] = 'protein_coding'
        # add in code to extract pacbio IDs here and add them to transcript biotype map
    consensus_dict = get_transcript_dict(args.consensus)
    consensus_gene_dict = get_gene_dict(consensus_dict)
    complete_reference_ids = extract_reference_completeness(args.refGff3)
    complete_ids = set.union(*[complete_reference_ids, complete_tm_ids, complete_aug_ids, complete_cgp_ids])
    basic = []
    for tx_list in consensus_gene_dict.itervalues():
        basic_coding = get_basic_coding(tx_list, complete_ids, transcript_biotype_map)
        if len(basic_coding) > 0:
            basic.extend(basic_coding)
            continue
        basic_noncoding = get_basic_noncoding(tx_list, complete_ids, transcript_biotype_map)
        if len(basic_noncoding) > 0:
            basic.extend(basic_noncoding)
            continue
        basic_pseudo = get_basic_pseudogenes(tx_list, transcript_biotype_map)
        if len(basic_pseudo) > 0:
            basic.extend(basic_pseudo)
            continue
        basic.append(sorted(tx_list, key=len, reverse=True)[0])
    ensureFileDir(args.outBasic)
    with open(args.outBasic, 'w') as outf:
        for tx in basic:
            outf.write('\t'.join(map(str, tx.get_gene_pred())) + '\n')


if __name__ == '__main__':
    main()
