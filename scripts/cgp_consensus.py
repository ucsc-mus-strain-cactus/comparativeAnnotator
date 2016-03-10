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
from collections import defaultdict, OrderedDict
from pycbio.sys.sqliteOps import open_database, get_multi_index_query_dict, get_query_dict
from pycbio.sys.fileOps import ensureFileDir
from comparativeAnnotator.database_queries import get_gene_transcript_map, get_transcript_gene_map
from pycbio.bio.transcripts import get_transcript_dict

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
        interval_map[x.chromosome].append(x.get_interval())
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
    gps = [consensus_dict[x] for x in ens_ids if x in consensus_dict]
    assert len(gps) != 0, cgp_tx.name
    # ignore the case where the CGP transcript was assigned transcripts on another chromosome
    full_intervals = {x for x in build_full_gene_intervals(gps)}
    for support, intron_interval in zip(*[intron_vector, cgp_tx.intron_intervals]):
        if not any([intron_interval.overlap(x) for x in full_intervals]) and support == 1:
            return True
    return False


def determine_if_new_introns(cgp_tx, ens_ids, consensus_dict, intron_vector):
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


def find_missing_transcripts(cgp_dict, consensus_genes, intron_dict, final_consensus, metrics, consensus_dict,
                             gene_transcript_map, support_cutoff=80.0):
    """
    If a CGP transcript is associated with genes that are all missing from the consensus, include it if it has at least
    support_cutoff supported introns. Otherwise, remove it. 
    Also, check whether this transcript does not at all overlap the consensus transcript(s) for the gene(s) it is
    assigned. If this is true, also re-include it. This is generally the result of transMap mapping two places.
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
            continue
        # does this transcript exist on a different chromosome from the consensus picked?
        gps = []
        for gene_id in cgp_genes:
            # gene may not be in gene_transcript map if it is listed as protein_coding but none of its transcripts are
            if gene_id not in gene_transcript_map:
                continue
            gps.extend([consensus_dict[x] for x in gene_transcript_map[gene_id] if x in consensus_dict])  
        full_intervals = build_full_gene_intervals(gps)
        cgp_interval = cgp_tx.get_interval()
        if not any([cgp_interval.overlap(x) for x in full_intervals]):
            percent_support = 100.0 * sum(intron_dict[cgp_id]) / len(intron_dict[cgp_id])
            if percent_support >= support_cutoff:
                final_consensus[cgp_id] = cgp_tx
                jg_genes.add(cgp_id.split(".")[0])
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
            cgp_tx, gene_id = replace_map[consensus_id]
            cgp_tx.id = cgp_tx.name
            cgp_tx.name = consensus_id
            cgp_tx.name2 = gene_id
            final_consensus[consensus_id] = cgp_tx
        else:
            final_consensus[consensus_id] = consensus_tx
    for cgp_tx in new_isoforms:
        final_consensus[cgp_tx.name] = cgp_tx


def update_transcripts(cgp_dict, consensus_dict, transcript_gene_map, intron_dict, final_consensus, metrics,
                       cgp_stats_dict, consensus_stats_dict):
    """
    Main transcript replacement/inclusion algorithm.
    For every cgp transcript, determine if it should replace one or more consensus transcripts.
    If it should not, then determine if it should be kept because it adds new splice junctions.
    """
    replace_map = {}  # will store a mapping between consensus IDs and the CGP transcripts that will replace them
    new_isoforms = []  # will store cgp transcripts which represent new potential isoforms of a gene
    join_genes = {"Unsupported": 0, "Supported": 0}  # will count the number of unsupported joins
    for cgp_id, cgp_tx in cgp_dict.iteritems():
        cgp_stats = cgp_stats_dict[cgp_id]
        ens_ids = {x for x in cgp_stats.keys() if x in consensus_stats_dict and consensus_dict[x].chromosome ==
                   cgp_tx.chromosome}  # we filter out join-gene cases that are on alternative chromosomes
        # as a result of this filtering, we cannot rely on name2 to show if this is a join-genes case
        gene_ids = {transcript_gene_map[x] for x in ens_ids}
        consensus_stats = {x: consensus_stats_dict[x] for x in ens_ids}
        to_replace_ids = determine_if_better(cgp_stats, consensus_stats)
        intron_vector = intron_dict[cgp_id]
        if len(to_replace_ids) > 0:
            for to_replace_id in to_replace_ids:
                gene_id = transcript_gene_map[to_replace_id]
                replace_map[to_replace_id] = [cgp_tx, gene_id]
        elif determine_if_new_introns(cgp_tx, ens_ids, consensus_dict, intron_vector) is True:
            # make sure this isn't joining two genes in an unsupported way
            if len(gene_ids) == 1:
                new_isoforms.append(cgp_tx)
            elif determine_if_split_gene_is_supported(consensus_dict, cgp_tx, ens_ids, intron_vector):
                new_isoforms.append(cgp_tx)
                join_genes["Supported"] += 1
            else:
                join_genes["Unsupported"] += 1
    # calculate some metrics for plots once all genomes are analyzed
    collapse_rate = len(set(zip(*replace_map.values())[0]))
    metrics["CgpReplace"] = {"CgpReplaceRate": len(replace_map), "CgpCollapseRate": collapse_rate}
    metrics["NewIsoforms"] = len(new_isoforms)
    metrics["JoinGeneSupported"] = join_genes
    build_final_consensus(consensus_dict, replace_map, new_isoforms, final_consensus)


def evaluate_cgp_consensus(consensus_dict, metrics):
    tx_names = OrderedDict((("transMap",  set()), ("AugustusTMR", set()), ("CGP", set())))
    gene_names = OrderedDict((("GENCODE", set()), ("CGP", set())))
    for cgp_tx in consensus_dict.itervalues():
        if "jg" in cgp_tx.name:
            tx_names["CGP"].add(cgp_tx.name)
            gene_names["CGP"].update([x for x in cgp_tx.name2.split(",") if 'jg' in x])
            gene_names["GENCODE"].update([x for x in cgp_tx.name2.split(",") if not 'jg' in x])
        elif "aug" in cgp_tx.id:
            tx_names["AugustusTMR"].add(cgp_tx.id)
            gene_names["GENCODE"].add(cgp_tx.name2)
        else:
            tx_names["transMap"].add(cgp_tx.id)
            gene_names["GENCODE"].add(cgp_tx.name2)
    transcript_stats = OrderedDict([[x, len(y)] for x, y in tx_names.iteritems()])
    gene_stats = OrderedDict([[x, len(y)] for x, y in gene_names.iteritems()])
    metrics["ConsensusStats"] = {"Transcript": transcript_stats, "Gene": gene_stats}


def cgp_consensus(args):
    gene_transcript_map = get_gene_transcript_map(args.ref_genome, args.comp_db, biotype="protein_coding")
    transcript_gene_map = get_transcript_gene_map(args.ref_genome, args.comp_db, biotype="protein_coding")
    # open CGP database -- we don't need comparativeAnnotator databases anymore
    con, cur = open_database(args.cgp_db)
    # load both consensus and CGP into dictionaries
    consensus_dict = get_transcript_dict(args.consensus_gp)
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
    final_consensus = {}
    metrics = {}
    # save all CGP transcripts which have no associated genes
    find_new_transcripts(cgp_dict, final_consensus, metrics)
    # save all CGP transcripts whose associated genes are not in the consensus
    consensus_genes = {x.name2 for x in consensus_dict.itervalues()}
    find_missing_transcripts(cgp_dict, consensus_genes, intron_dict, final_consensus, metrics, consensus_dict,
                             gene_transcript_map)
    # remove all such transcripts from the cgp dict before we evaluate for updating
    cgp_dict = {x: y for x, y in cgp_dict.iteritems() if x not in final_consensus}
    update_transcripts(cgp_dict, consensus_dict, transcript_gene_map, intron_dict, final_consensus, metrics,
                       cgp_stats_dict, consensus_stats_dict)
    evaluate_cgp_consensus(final_consensus, metrics)
    # write results out to disk
    ensureFileDir(args.metrics_dir)
    with open(os.path.join(args.metrics_dir, args.genome + ".metrics.pickle"), "w") as outf:
        pickle.dump(metrics, outf)
    s = sorted(final_consensus.itervalues(), key=lambda tx: (tx.chromosome, tx.start))
    return ['\t'.join(map(str, x.get_gene_pred())) + '\n' for x in s]
