"""
Takes a comparative Augustus transcript set with the name2 field set to the best gencode gene ID.
Produces a new consensus from this transcript set, following these rules:

1) Align every CGP transcript to the consensus transcript's reference transcript. If one or more CGP transcript 
fulfills the coverage/identity heuristics above, replace the consensus transcript with the best CGP by % identity.
2) If there are any remaining CGP transcripts that have not been moved into the consensus set, we look for RNAseq 
supported splice junctions not present in any of the transcripts for this gene. If so, we include it. 
3) If any augustus CGP transcripts do not overlap any existing transcripts, include them

"""
import argparse
import os
import cPickle as pickle
from collections import OrderedDict, defaultdict
from pycbio.sys.sqliteOps import open_database, get_multi_index_query_dict, get_query_dict
from pycbio.sys.fileOps import ensureDir
from comparativeAnnotator.database_queries import get_transcript_gene_map, get_transcript_biotype_map
from pycbio.bio.transcripts import get_transcript_dict, Transcript

__author__ = "Ian Fiddes"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cgpDb", default="cgp_cds_metrics.db")
    parser.add_argument("--cgpGp", required=True)
    parser.add_argument("--consensusGp", required=True)
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--intronBitsPath", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--refGenome", required=True)
    parser.add_argument("--outGp", required=True)
    parser.add_argument("--metricsOutDir", required=True)
    return parser.parse_args()


def resolve_multiple_parents(cgp_stats_dict, cgp_dict, transcript_biotype_map, transcript_gene_map,
                             cov_cutoff=0.35, ident_cutoff=0.90, cov_weight=0.25, ident_weight=0.75):
    """
    Resolves CGP transcripts which were assigned to more than one gene by picking a gene based on which transcripts
    had the highest score.
    """
    for cgp_id, cgp_tx in cgp_dict.iteritems():
        if ',' not in cgp_tx.name2:
            continue
        cgp_stats = cgp_stats_dict[cgp_id]
        scores = {}
        for ens_id, (cgp_cov, cgp_ident) in cgp_stats.iteritems():
            if cgp_cov <= cov_cutoff or cgp_ident <= ident_cutoff:
                continue
            elif transcript_biotype_map[ens_id] != 'protein_coding':
                continue
            score = cgp_cov * cov_weight + cgp_ident * ident_weight
            scores[ens_id] = score
        if len(scores) == 0:
            cgp_tx.name2 = cgp_id.split('.')[0]  # remove association with any gene
        else:
            best_tx = sorted(scores.iteritems(), key=lambda x: -x[1])[0][0]
            cgp_tx.name2 = transcript_gene_map[best_tx]


def load_intron_bits(intron_bits_path):
    """
    Load the intron bit vector files into a dictionary, properly handling cases where there are no introns
    """
    intron_dict = {}
    for line in open(intron_bits_path):
        l = line.split()
        if len(l) == 1:
            intron_dict[l[0]] = []
        else:
            intron_dict[l[0]] = map(int, l[1].split(","))
    return intron_dict


def build_splice_junction_set(gps):
    """
    Given an iterable of GenePredTranscript objects, returns a set of all splice junction intervals
    as ChromosomeIntervals
    """
    sjs = set()
    for gp in gps:
        for intron_interval in gp.intron_intervals:
            sjs.add(intron_interval)
    return sjs


def filter_cgp_splice_junctions(gp, intron_vector):
    """
    Returns a set of ChromosomeInterval objects that are filtered based on the intron vector provided
    """
    sjs = set()
    assert len(intron_vector) == len(gp.intron_intervals)
    for support, intron_interval in zip(*[intron_vector, gp.intron_intervals]):
        if support == 1:
            sjs.add(intron_interval)
    return sjs


def build_full_gene_intervals(gps):
    """
    Given a iterable of GenePredTranscripts, return ChromosomeIntervals that encapsulate the maximum boundaries of
    these intervals. This can be more than one if the gene ends up mapping to different chromosomes.
    """
    interval_map = defaultdict(list)
    for x in gps:
        interval_map[x.chromosome].append(x.get_interval())
    r = set()
    for intervals in interval_map.itervalues():
        r.add(reduce(lambda x, y: x.hull(y), intervals))
    return r


def determine_if_new_introns(cgp_tx, ens_ids, tmr_consensus_dict, intron_vector):
    """
    Use intron bit information to build a set of CGP introns, and compare this set to all consensus introns
    for a given set of genes.
    """
    gps = [tmr_consensus_dict[x] for x in ens_ids if x in tmr_consensus_dict]
    ens_splice_junctions = build_splice_junction_set(gps)
    cgp_splice_junctions = filter_cgp_splice_junctions(cgp_tx, intron_vector)
    if len(cgp_splice_junctions - ens_splice_junctions) > 0:
        return True
    return False


def determine_if_better(cgp_stats, consensus_stats, cov_weight=0.25, ident_weight=0.75):
    """
    Determines if this CGP transcript is better than any of the consensus transcripts it may come from
    """
    ens_ids = []
    for ens_id, (consensus_cov, consensus_ident) in consensus_stats.iteritems():
        cgp_cov, cgp_ident = cgp_stats[ens_id]
        if (cov_weight * cgp_cov + ident_weight * cgp_ident) > (cov_weight * consensus_cov + ident_weight * consensus_ident):
            ens_ids.append(ens_id)
    return ens_ids


def find_new_transcripts(cgp_dict, consensus, metrics):
    """
    Include any transcripts which were not assigned any gene IDs
    """
    jg_genes = set()
    for cgp_id, cgp_tx in cgp_dict.iteritems():
        if 'jg' in cgp_tx.name2:
            consensus[cgp_id] = cgp_tx
            jg_genes.add(cgp_id.split(".")[0])
    metrics["CgpAdditions"] = {"CgpNewGenes": len(jg_genes), "CgpNewTranscripts": len(consensus)}


def find_missing_transcripts(cgp_dict, consensus_genes, intron_dict, consensus, metrics, support_cutoff=80.0):
    """
    If a CGP transcript is associated with a gene missing from the consensus, include it.
    """
    jg_genes = set()
    to_remove = set()
    for cgp_id, cgp_tx in cgp_dict.iteritems():
        if 'jg' in cgp_tx.name2:
            continue
        cgp_genes = set(cgp_tx.name2.split(","))
        if len(cgp_genes & consensus_genes) == 0:
            percent_support = 100.0 * sum(intron_dict[cgp_id]) / len(intron_dict[cgp_id])
            if percent_support >= support_cutoff:
                consensus[cgp_id] = cgp_tx
                jg_genes.add(cgp_id.split(".")[0])
            # we want to exclude transcripts without support from further consensus finding
            else:
                to_remove.add(cgp_id)
            continue
    metrics["CgpAddMissing"] = {"CgpMissingGenes": len(jg_genes), "CgpMissingTranscripts": len(consensus)}
    # hacking this in here
    for cgp_id in to_remove:
        del cgp_dict[cgp_id]


def build_consensus(tmr_consensus_dict, replace_map, new_isoforms, consensus):
    """
    Builds the final consensus gene set given the replace map as well as the new isoforms. Deduplicates based on id.
    """
    seen_ids = set()  # used to deduplicate the results
    for consensus_id, consensus_tx in tmr_consensus_dict.iteritems():
        if consensus_id in replace_map and consensus_id not in seen_ids:
            seen_ids.add(consensus_id)
            cgp_tx, gene_id = replace_map[consensus_id]
            cgp_tx.id = cgp_tx.name
            cgp_tx.name = consensus_id
            cgp_tx.name2 = gene_id
            consensus[consensus_id] = cgp_tx
        else:
            consensus[consensus_id] = consensus_tx
    for cgp_tx in new_isoforms:
        consensus[cgp_tx.name] = cgp_tx


def update_transcripts(cgp_dict, tmr_consensus_dict, transcript_gene_map, intron_dict, consensus, metrics,
                       cgp_stats_dict, consensus_stats_dict, transcript_biotype_map):
    """
    Main transcript replacement/inclusion algorithm.
    For every cgp transcript, determine if it should replace one or more consensus transcripts.
    If it should not, then determine if it should be kept because it adds new splice junctions.
    """
    replace_map = {}  # will store a mapping between consensus IDs and the CGP transcripts that will replace them
    new_isoforms = []  # will store cgp transcripts which represent new potential isoforms of a gene
    for cgp_id, cgp_tx in cgp_dict.iteritems():
        cgp_stats = cgp_stats_dict[cgp_id]
        # we don't want to try and replace non-coding transcripts with CGP predictions
        ens_ids = {x for x in cgp_stats.keys() if transcript_biotype_map[x] == 'protein_coding'}
        consensus_stats = {x: consensus_stats_dict[x] for x in ens_ids if x in consensus_stats_dict}
        to_replace_ids = determine_if_better(cgp_stats, consensus_stats)
        intron_vector = intron_dict[cgp_id]
        if len(to_replace_ids) > 0:
            for to_replace_id in to_replace_ids:
                gene_id = transcript_gene_map[to_replace_id]
                replace_map[to_replace_id] = [cgp_tx, gene_id]
        elif determine_if_new_introns(cgp_tx, ens_ids, tmr_consensus_dict, intron_vector) is True:
            new_isoforms.append(cgp_tx)
    # calculate some metrics for plots once all genomes are analyzed
    collapse_rate = len(set(zip(*replace_map.values())[0]))
    metrics["CgpReplace"] = {"CgpReplaceRate": len(replace_map), "CgpCollapseRate": collapse_rate}
    metrics["NewIsoforms"] = len(new_isoforms)
    build_consensus(tmr_consensus_dict, replace_map, new_isoforms, consensus)


def deduplicate_cgp_consensus(cgp_consensus, metrics):
    """
    In the process of consensus building, we may find that we have ended up with more than one transcript for a gene
    that are actually identical. Remove these, picking the non-CGP gene in these cases.
    """
    deduplicated_consensus = []
    duplicates = defaultdict(lambda: defaultdict(list))
    for tx in cgp_consensus.itervalues():
        if tx.cds_size <= 25:
            deduplicated_consensus.append(tx)
        else:
            # generate a BED of only CDS coordinates
            tx_cds_bed = tx.get_bed(start_offset=tx.thick_start, stop_offset=tx.thick_stop)
            # create a new transcript object out of this
            tx_cds = Transcript(tx_cds_bed)
            duplicates[tx.name2][frozenset(tx_cds.exon_intervals)].append(tx)
    dup_count = 0
    for gene in duplicates:
        for cds_interval, tx_list in duplicates[gene].iteritems():
            cgp_txs = [tx for tx in tx_list if 'jg' in tx.id]
            if len(cgp_txs) > 0:
                dup_count += 1
            else:
                deduplicated_consensus.extend(tx_list)
    metrics['CgpEqualsTMR'] = dup_count
    return deduplicated_consensus


def evaluate_cgp_consensus(consensus, metrics):
    tx_names = OrderedDict((("transMap",  set()), ("AugustusTMR", set()), ("CGP", set())))
    gene_names = OrderedDict((("Gencode", set()), ("CGP", set())))
    for cgp_tx in consensus:
        if "jg" in cgp_tx.name:
            tx_names["CGP"].add(cgp_tx.name)
            gene_names["CGP"].update([x for x in cgp_tx.name2.split(",") if 'jg' in x])
            gene_names["Gencode"].update([x for x in cgp_tx.name2.split(",") if 'jg' not in x])
        elif "aug" in cgp_tx.id:
            tx_names["AugustusTMR"].add(cgp_tx.id)
            gene_names["Gencode"].add(cgp_tx.name2)
        else:
            tx_names["transMap"].add(cgp_tx.id)
            gene_names["Gencode"].add(cgp_tx.name2)
    transcript_stats = OrderedDict([[x, len(y)] for x, y in tx_names.iteritems()])
    gene_stats = OrderedDict()
    gene_stats['Gencode'] = len(gene_names['Gencode'] - gene_names['CGP'])
    gene_stats['Both'] = len(gene_names['Gencode'] & gene_names['CGP'])
    gene_stats['CGP'] = len(gene_names['CGP'] - gene_names['Gencode'])
    metrics["ConsensusStats"] = {"Transcript": transcript_stats, "Gene": gene_stats}


def cgp_consensus(args):
    transcript_biotype_map = get_transcript_biotype_map(args.ref_genome, args.comp_db)
    transcript_gene_map = get_transcript_gene_map(args.ref_genome, args.comp_db, biotype="protein_coding")
    con, cur = open_database(args.cgp_db)
    # load both consensus and CGP into dictionaries
    tmr_consensus_dict = get_transcript_dict(args.consensus_gp)
    cgp_dict = get_transcript_dict(args.cgp_gp)
    # load the BLAT results from the sqlite database
    cgp_stats_query = "SELECT CgpId,EnsId,AlignmentCoverage,AlignmentIdentity FROM '{}_cgp'".format(args.genome)
    cgp_stats_dict = get_multi_index_query_dict(cur, cgp_stats_query, num_indices=2)
    consensus_stats_query = ("SELECT EnsId,AlignmentCoverage,AlignmentIdentity FROM "
                             "'{}_consensus'".format(args.genome))
    consensus_stats_dict = get_query_dict(cur, consensus_stats_query)
    # load the intron bits
    intron_dict = load_intron_bits(args.cgp_intron_bits)
    # final dictionaries
    cgp_consensus = {}
    metrics = {}
    # resolve transcripts assigned to multiple genes
    resolve_multiple_parents(cgp_stats_dict, cgp_dict, transcript_biotype_map, transcript_gene_map)
    # save all CGP transcripts which have no associated genes
    find_new_transcripts(cgp_dict, cgp_consensus, metrics)
    # save all CGP transcripts whose associated genes are not in the consensus
    consensus_genes = {x.name2 for x in tmr_consensus_dict.itervalues()}
    find_missing_transcripts(cgp_dict, consensus_genes, intron_dict, cgp_consensus, metrics)
    # remove all such transcripts from the cgp dict before we evaluate for updating
    cgp_dict = {x: y for x, y in cgp_dict.iteritems() if x not in cgp_consensus}
    update_transcripts(cgp_dict, tmr_consensus_dict, transcript_gene_map, intron_dict, cgp_consensus, metrics,
                       cgp_stats_dict, consensus_stats_dict, transcript_biotype_map)
    deduplicated_consensus = deduplicate_cgp_consensus(cgp_consensus, metrics)
    evaluate_cgp_consensus(deduplicated_consensus, metrics)
    # write results out to disk
    ensureDir(args.metrics_dir)
    with open(os.path.join(args.metrics_dir, args.genome + ".metrics.pickle"), "w") as outf:
        pickle.dump(metrics, outf)
    s = sorted(deduplicated_consensus, key=lambda tx: (tx.chromosome, tx.start))
    return ['\t'.join(map(str, x.get_gene_pred())) + '\n' for x in s]
