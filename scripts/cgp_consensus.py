"""
Takes a comparative Augustus transcript set with the name2 field set to the best gencode gene ID.
Produces a new consensus from this transcript set, following these rules:

1) Align every CGP transcript to the consensus transcript's reference transcript. If one or more CGP transcript 
fufills the coverage/identity heuristics above, replace the consensus transcript with the best CGP by % identity.
2) If there are any remaining CGP transcripts that have not been moved into the consensus set, we look for RNAseq 
supported splice junctions not present in any of the transcripts for this gene. If so, we include it. 
3) If any augustus CGP transcripts do not overlap any existing transcripts, include them

"""
import argparse
import os
import cPickle as pickle
from collections import defaultdict
import lib.sql_lib as sql_lib
import lib.psl_lib as psl_lib
import lib.seq_lib as seq_lib

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
        interval_map[x.chromosome].append(seq_lib.ChromosomeInterval(x.chromosome, x.start, x.stop, x.strand))
    r = set()
    for intervals in interval_map.itervalues():
        r.add(reduce(lambda x, y: x.hull(y), intervals))
    return r


def determine_if_split_gene_is_supported(consensus_dict, cgp_tx, ens_ids, intron_vector):
    """
    If a CGP gene is assigned more than one gene, determine if the spanning intron blocks are supported by RNAseq.
    Returns True if this joined gene is supported
    A canonical example:
    -----------------CGP----------------
    -----ENSG1-----         -----ENSG2-----
    ----ENST1A-----         ----ENST2A-----
    --ENST1B--                  --ENST2B--
    """
    cgp_genes = cgp_tx.name2.split(",")
    full_intervals = set()
    for cgp_gene in cgp_genes:
        txs = [consensus_dict[x] for x in ens_ids if x in consensus_dict]
        assert len(txs) != 0, cgp_tx.name
        full_intervals.update(build_full_gene_intervals(txs))
    for support, intron_interval in zip(*[intron_vector, cgp_tx.intron_intervals]):
        if not any([intron_interval.overlap(x) for x in full_intervals]) and support == 1:
            return True
    return False


def determine_if_new_introns(cgp_id, cgp_tx, ens_ids, consensus_dict, intron_vector):
    """
    Use intron bit information to build a set of CGP introns, and compare this set to all consensus introns
    for a given set of genes.
    """
    gps = [consensus_dict[x] for x in ens_ids if x in consensus_dict]
    ens_splice_junctions = build_splice_junction_set(gps)
    cgp_splice_junctions = filter_cgp_splice_junctions(cgp_tx, intron_vector)
    if len(cgp_splice_junctions - ens_splice_junctions) > 0:
        return True
    return False


def determine_if_better(cgp_stats, consensus_stats):
    """
    Determines if this CGP transcript is better than any of the consensus transcripts it may come from
    """
    ens_ids = []
    for ens_id, (consensus_cov, consensus_ident) in consensus_stats.iteritems():
        cgp_cov, cgp_ident = cgp_stats[ens_id]
        if ((cgp_ident > consensus_ident and cgp_cov >= consensus_cov) or 
                (cgp_cov > consensus_cov and cgp_ident >= consensus_ident)):
            ens_ids.append(ens_id)
    return ens_ids


def find_new_transcripts(cgp_dict, final_consensus, metrics):
    """
    Include any transcripts which were not assigned any gene IDs
    """
    jg_genes = set()
    for cgp_id, cgp_tx in cgp_dict.iteritems():
        if 'jg' in cgp_tx.name2:
            final_consensus[cgp_id] = cgp_tx
            jg_genes.add(cgp_id.split(".")[0])
    metrics["CgpAdditions"] = {"CgpNewGenes": len(jg_genes), "CgpNewTranscripts": len(final_consensus)}


def find_missing_transcripts(cgp_dict, consensus_genes, intron_dict, final_consensus, metrics, support_cutoff=80.0):
    """
    If a CGP transcript is associated with genes that are all missing from the consensus, include it if it has at least
    support_cutoff supported introns. Otherwise, remove it.
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
                final_consensus[cgp_id] = cgp_tx
                jg_genes.add(cgp_id.split(".")[0])
            # we want to exclude transcripts without support from further consensus finding
            else:
                to_remove.add(cgp_id)
    metrics["CgpAddMissing"] = {"CgpMissingGenes": len(jg_genes), "CgpMissingTranscripts": len(final_consensus)}
    # hacking this in here
    for cgp_id in to_remove:
        del cgp_dict[cgp_id]


def build_final_consensus(consensus_dict, replace_map, new_isoforms, final_consensus):
    """
    Builds the final consensus gene set given the replace map as well as the new isoforms. Deduplicates.
    """
    seen_ids = set()  # used to deduplicate the results
    for consensus_id, consensus_tx in consensus_dict.iteritems():
        if consensus_id in replace_map and consensus_id not in seen_ids:
            seen_ids.add(consensus_id)
            cgp_tx = replace_map[consensus_id]
            cgp_tx.id = cgp_tx.name
            cgp_tx.name = consensus_id
            final_consensus[consensus_id] = cgp_tx
        else:
            final_consensus[consensus_id] = consensus_tx
    for cgp_tx in new_isoforms:
        final_consensus[cgp_tx.name] = cgp_tx


def update_transcripts(cgp_dict, consensus_dict, genome, gene_transcript_map, intron_dict, final_consensus, metrics,
                       cgp_stats_dict, consensus_stats_dict):
    """
    Main transcript replacement/inclusion algorithm.
    For every cgp transcript, determine if it should replace one or more consensus transcripts.
    If it should not, then determine if it should be kept because it adds new splice junctions.
    """
    replace_map = {}  # will store a mapping between consensus IDs and the CGP IDs that will replace them
    new_isoforms = []  # will store cgp IDs which represent new potential isoforms of a gene
    join_genes = {"Unsupported": 0, "Supported": 0}  # will count the number of unsupported joins
    for cgp_id, cgp_tx in cgp_dict.iteritems():
        cgp_stats = cgp_stats_dict[cgp_id]
        ens_ids = {x for x in cgp_stats.keys() if x in consensus_stats_dict}
        consensus_stats = {x: consensus_stats_dict[x] for x in ens_ids}
        to_replace_ids = determine_if_better(cgp_stats, consensus_stats)
        intron_vector = intron_dict[cgp_id]
        if len(to_replace_ids) > 0:
            for to_replace_id in to_replace_ids:
                replace_map[to_replace_id] = cgp_tx
        elif determine_if_new_introns(cgp_id, cgp_tx, ens_ids, consensus_dict, intron_vector) is True:
            # make sure this isn't joining two genes in an unsupported way
            if len(cgp_tx.name2.split(",")) == 1:
                new_isoforms.append(cgp_tx)
            elif determine_if_split_gene_is_supported(consensus_dict, cgp_tx, ens_ids, intron_vector):
                new_isoforms.append(cgp_tx)
                join_genes["Supported"] += 1
            else:
                join_genes["Unsupported"] += 1
    # calculate some metrics for plots once all genomes are analyzed
    metrics["CgpReplace"] = {"CgpReplaceRate": len(replace_map), "CgpCollapseRate": len(set(replace_map.itervalues()))}
    metrics["NewIsoforms"] = len(new_isoforms)
    metrics["JoinGeneSupported"] = join_genes
    build_final_consensus(consensus_dict, replace_map, new_isoforms, final_consensus)


def main():
    args = parse_args()
    # attach regular comparativeAnnotator reference databases in order to build gene-transcript map
    con, cur = sql_lib.attach_databases(args.compAnnPath, mode="reference")
    gene_transcript_map = sql_lib.get_gene_transcript_map(cur, args.refGenome, biotype="protein_coding")
    # open CGP database -- we don't need comparativeAnnotator databases anymore
    cgp_db = os.path.join(args.compAnnPath, args.cgpDb)
    con, cur = sql_lib.open_database(cgp_db)
    # load both consensus and CGP into dictionaries
    consensus_dict = seq_lib.get_transcript_dict(args.consensusGp)
    cgp_dict = seq_lib.get_transcript_dict(args.cgpGp)
    # load the BLAT results from the sqlite database
    cgp_stats_query = "SELECT CgpId,EnsId,AlignmentCoverage,AlignmentIdentity FROM '{}_cgp'".format(args.genome)
    cgp_stats_dict = sql_lib.get_multi_index_query_dict(cur, cgp_stats_query, num_indices=2)
    consensus_stats_query = ("SELECT EnsId,AlignmentCoverage,AlignmentIdentity FROM "
                             "'{}_consensus'".format(args.genome))
    consensus_stats_dict = sql_lib.get_query_dict(cur, consensus_stats_query)
    # load the intron bits
    intron_dict = load_intron_bits(args.intronBitsPath)
    # final dictionaries
    final_consensus = {}
    metrics = {}
    # save all CGP transcripts which have no associated genes
    find_new_transcripts(cgp_dict, final_consensus, metrics)
    # save all CGP transcripts whose associated genes are not in the consensus
    consensus_genes = {x.name2 for x in consensus_dict.itervalues()}
    find_missing_transcripts(cgp_dict, consensus_genes, intron_dict, final_consensus, metrics)
    # remove all such transcripts from the cgp dict before we evaluate for updating
    cgp_dict = {x: y for x, y in cgp_dict.iteritems() if x not in final_consensus}
    update_transcripts(cgp_dict, consensus_dict, args.genome, gene_transcript_map, intron_dict, final_consensus, 
                       metrics, cgp_stats_dict, consensus_stats_dict)
    # write results out to disk
    with open(os.path.join(args.metricsOutDir, args.genome + ".metrics.pickle"), "w") as outf:
        pickle.dump(metrics, outf)
    with open(args.outGp, "w") as outf:
        for tx_id, tx in sorted(final_consensus.iteritems(), key=lambda x: [x[1].chromosome, x[1].start]):
            outf.write("\t".join(map(str, tx.get_gene_pred())) + "\n")


if __name__ == "__main__":
    main()
