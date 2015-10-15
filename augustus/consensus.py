import argparse
import re
import os
from collections import defaultdict, OrderedDict, Counter
import lib.sql_lib as sql_lib
import lib.psl_lib as psl_lib
from lib.general_lib import mkdir_p, merge_dicts
import etc.config


__author__ = "Ian Fiddes"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--refGenome", required=True)
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--augGp", required=True)
    parser.add_argument("--tmGp", required=True)
    return parser.parse_args()


def load_gps(gp_paths):
    """
    Get a dictionary mapping all gene IDs from a genePred into its entire record. If the gene IDs are not unique
    this function will not work like you want it to.
    """
    return {l.split()[0]: l for p in gp_paths for l in open(p)}


def merge_stats(cur, genome):
    """
    Adapter function to return the combination of the stats dicts produced for augustus and transMap
    """
    tm_stats = sql_lib.get_stats(cur, genome, mode="transMap")
    aug_stats = sql_lib.get_stats(cur, genome, mode="augustus")
    return merge_dicts([tm_stats, aug_stats])


def find_best_aln(stats, sig_fig=4):
    """
    Takes a list of stats for transcripts and finds all transcripts which have the highest percent identity
    """
    s = sorted(stats, key=lambda x: -x[1])
    best_ident = round(s[0][1], sig_fig)
    return {name for name, aln_id, aln_cov in s if round(aln_id, sig_fig) == best_ident}


def build_data_dict(id_names, id_list, transcript_gene_map):
    """
    Builds a dictionary mapping gene_id -> transcript_ids -> aln_ids in id_names bins (as an OrderedDict)
    """
    data_dict = defaultdict(dict)
    for ids, n in zip(*[id_list, id_names]):
        for aln_id in ids:
            ens_id = psl_lib.strip_alignment_numbers(aln_id)
            gene_id = transcript_gene_map[ens_id]
            if ens_id not in data_dict[gene_id]:
                data_dict[gene_id][ens_id] = OrderedDict((x, []) for x in id_names)
            data_dict[gene_id][ens_id][n].append(aln_id)
    return data_dict


def find_best_alns(stats, ids):
    """
    Takes the list of transcript Ids and finds the best alignment(s)
    """
    s = []
    for aln_id in ids:
        ident, cov = stats[aln_id]
        s.append([aln_id, ident, cov])
    s = sorted(s, key=lambda x: -x[1])
    best_val = round(s[0][1], 6)
    return [x[0] for x in s if round(x[1], 6) == best_val]


def evaluate_ids(fail_ids, good_specific_ids, pass_ids, aug_ids, stats):
    """
    For a given ensembl ID, we have augustus/transMap ids in 4 categories. Based on the hierarchy Pass>Good>Fail,
    return the best transcript in the highest category with a transMap transcript. 
    """
    if len(pass_ids) > 0:
        best_alns = find_best_alns(stats, pass_ids + aug_ids)
        return best_alns, "pass"
    elif len(good_specific_ids) > 0:
        best_alns = find_best_alns(stats, good_specific_ids + aug_ids)
        return best_alns, "good"
    elif len(fail_ids) > 0:
        best_alns = find_best_alns(stats, fail_ids + aug_ids)
        return best_alns, "fail"
    else:
        return None, "no_good_aln"


def is_tie(best_alns):
    """
    If we have more than one best transcript, is at least one from transMap and one from Augustus?
    """
    seen = set()
    for aln_id in best_alns:
        ens_id = psl_lib.remove_augustus_alignment_number(aln_id)
        if ens_id in seen:
            return True
        else:
            seen.add(ens_id)
    return False


def find_best_transcripts(data_dict, stats):
    """
    For all of the transcripts categorized in data_dict, evaluate them and bin them.
    """
    binned_transcripts = {}
    for gene_id in data_dict:
        binned_transcripts[gene_id] = {}
        for ens_id in data_dict[gene_id]:
            tx_recs = data_dict[gene_id][ens_id]
            fail_ids, good_specific_ids, pass_ids, aug_ids = tx_recs.values()
            best_alns, category = evaluate_ids(fail_ids, good_specific_ids, pass_ids, aug_ids, stats)
            if best_alns is None:
                binned_transcripts[gene_id][ens_id] = [best_alns, category, None]
            else:
                tie = is_tie(best_alns)
                binned_transcripts[gene_id][ens_id] = [best_alns[0], category, tie]
    return binned_transcripts


def find_best_for_gene(bins, stats, cov_cutoff=80.0, ident_cutoff=80.0):
    aln_ids = zip(*bins.itervalues())[0]
    keep_ids = []
    for aln_id in aln_ids:
        if aln_id is None:
            continue
        ident, cov = stats[aln_id]
        if cov >= cov_cutoff and ident >= ident_cutoff:
            keep_ids.append(aln_id)
    if len(keep_ids) > 0:
        return find_best_alns(stats, keep_ids)
    else:
        return None


def find_consensus(binned_transcripts, stats):
    """
    Takes the binned transcripts and picks a consensus. This works like this:
    For each gene, evalaute the transcripts. If the best alignment falls 
    """
    consensus = []
    for gene_id in binned_transcripts:
        gene_in_consensus = False
        for ens_id in binned_transcripts[gene_id]:
            best_id, category, is_tie = binned_transcripts[gene_id][ens_id]
            if category in ["pass", "good"]:
                consensus.append(best_id)
                gene_in_consensus = True
        if gene_in_consensus is False:
            # find the one best transcript for this gene based on identity after filtering for coverage
            best_for_gene = find_best_for_gene(binned_transcripts[gene_id], stats)
            if best_for_gene is not None:
                consensus.append(best_for_gene[0])
    return consensus


def consensus_by_biotype(cur, ref_genome, genome, biotype, transcript_gene_map, stats):
    fail_ids, good_specific_ids, pass_ids = sql_lib.get_fail_good_pass_ids(cur, ref_genome, genome, biotype)
    chr_y_ids = {psl_lib.strip_alignment_numbers(x) for x in sql_lib.get_ids_by_chromosome(cur, genome)}
    biotype_ids = sql_lib.get_biotype_ids(cur, ref_genome, biotype)
    # hacky way to avoid duplicating code in consensus finding - we will always have an aug_id set, it just may be empty
    if biotype == "protein_coding":
        aug_query = etc.config.augustusEval(genome)
        aug_ids = sql_lib.get_query_ids(cur, aug_query)
    else:
        aug_ids = set()
    id_names = ["fail_ids", "good_specific_ids", "pass_ids", "aug_ids"]
    id_list = [fail_ids, good_specific_ids, pass_ids, aug_ids]
    data_dict = build_data_dict(id_names, id_list, transcript_gene_map)
    binned_transcripts = find_best_transcripts(data_dict, stats)
    consensus = find_consensus(binned_transcripts, stats)
    return consensus


def fix_gene_pred(gp, transcript_gene_map):
    """
    These genePreds have a few problems. First, the alignment numbers must be removed. Second, we want to fix
    the name2 field to be the gene name. Third, we want to set the unique ID field. Finally, we want to sort the whole 
    thing by genomic coordinates.
    """
    gp = sorted([x.split("\t") for x in gp], key=lambda x: [x[1], x[3]])
    fixed = []
    for x in gp:
        x[10] = x[0]  # use unique Aug/TM ID as unique identifier
        tx_id = psl_lib.strip_alignment_numbers(x[0])
        x[0] = tx_id
        gene_id = transcript_gene_map[tx_id]
        x[11] = gene_id
        fixed.append(x)
    return ["\t".join(x) for x in fixed]


def write_gps(consensus, gps, consensus_base_path, biotype, transcript_gene_map):
    p = os.path.join(consensus_base_path, biotype + ".consensus_gene_set.gp")
    mkdir_p(os.path.dirname(p))
    try:
        gp_recs = [gps[aln_id] for aln_id in consensus]
    except:
        assert False, (consensus, biotype)
    fixed_gp_recs = fix_gene_pred(gp_recs, transcript_gene_map)
    with open(p, "w") as outf:
        for rec in fixed_gp_recs:
            outf.write(rec + "\n")


def main():
    args = parse_args()
    con, cur = sql_lib.attach_databases(args.compAnnPath, mode="augustus")
    transcript_gene_map = sql_lib.get_transcript_gene_map(cur, args.refGenome)
    biotypes = sql_lib.get_all_biotypes(cur, args.refGenome, gene_level=True)
    gps = load_gps([args.tmGp, args.augGp])  # load all Augustus and transMap transcripts into one big dict
    consensus_base_path = os.path.join(args.outDir, "gene_sets", args.genome)
    stats = merge_stats(cur, args.genome)
    for biotype in biotypes:
        consensus = consensus_by_biotype(cur, args.refGenome, args.genome, biotype, transcript_gene_map, stats)
        if len(consensus) > 0:  # some biotypes we may have nothing
            write_gps(consensus, gps, consensus_base_path, biotype, transcript_gene_map)



if __name__ == "__main__":
    main()