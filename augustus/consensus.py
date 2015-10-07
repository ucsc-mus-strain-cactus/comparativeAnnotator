import argparse
import re
import os
import itertools
from collections import defaultdict
from plotting.plot_functions import get_all_biotypes, get_gene_map, get_gene_biotype_map, gp_chrom_filter, get_all_ok, \
                                    load_gps, get_all_ids, get_reverse_name_map, strip_alignment_numbers, transmap_ok
import cPickle as pickle
import lib.sql_lib as sql_lib
from lib.general_lib import mkdir_p


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--binnedTranscriptPath", required=True)
    parser.add_argument("--attributePath", required=True)
    parser.add_argument("--augGp", required=True)
    parser.add_argument("--tmGp", required=True)
    parser.add_argument("--compGp", required=True)
    parser.add_argument("--basicGp", required=True)
    return parser.parse_args()


def merge_stats(cur, genome):
    """
    Adapter function to return the combination of the stats dicts produced for augustus and transMap
    """
    tm_stats = sql_lib.get_stats(cur, genome, category="transMap")
    aug_stats = sql_lib.get_stats(cur, genome, category="augustus")
    r = tm_stats.copy()
    r.update(aug_stats)
    return r


def find_ok_not_ok_candidates(stats_dict, ids, ok_ids, discard_cov_cutoff, filter_cov_cutoff):
    """
    Returns a list of candidate transcripts for OK and not OK. All transcripts must have discard_cov_cutoff
    coverage or they are discarded. In order to be a OK candidate, a transcript must be classifier OK and have
    coverage above filter_cov_cutoff.
    """
    ok_candidates = []
    not_ok_candidates = []
    discarded = []
    for aln_id in ids:
        if stats_dict[aln_id][2] < discard_cov_cutoff:
            discarded.append(stats_dict[aln_id][0])
        elif aln_id in ok_ids and stats_dict[aln_id][2] >= filter_cov_cutoff:
            ok_candidates.append(stats_dict[aln_id])
        else:
            not_ok_candidates.append(stats_dict[aln_id])
    return ok_candidates, not_ok_candidates, discarded


def find_best_aln(stats):
    """
    Takes the list of OK/notOK candidates from find_ok_not_ok_candidates and returns the best alignment(s).
    """
    s = sorted(stats, key=lambda x: -x[1])
    best_ident = round(s[0][1], 6)
    return [x[0] for x in s if round(x[1], 6) == best_ident]


def split_alternatives(stats, best_id):
    """
    Takes the list of OK/notOK candidates and splits them up by augustus/transMap, filtering out the winner
    """
    tm = []
    aug = []
    for aln_id, ident, cov in stats:
        if aln_id == best_id:
            continue
        if aln_id.startswith("aug"):
            aug.append(aln_id)
        else:
            tm.append(aln_id)
    return aug, tm


def analyze_candidates(candidates, aug_re=re.compile("I[1-2]-")):
    """
    Analyzes candidate transcripts, finding the winner and splitting the alternatives based on coming from augustus
    or transMap. If there are multiple equal winners from both transMap and Augustus, reports so.
    aug_re is used if Augustus was ran in two modes (I1 and I2) to prevent counting that as a tie with transMap.
    """
    winner_ids = find_best_aln(candidates)
    winner_id = winner_ids[0]
    if len({"-".join(aug_re.split(x[0])) for x in winner_ids}) > 1:
        is_tie = True
    else:
        is_tie = False
    aug_alts, tm_alts = split_alternatives(candidates, winner_id)
    return winner_id, aug_alts, tm_alts, is_tie


def bin_candidates(winner_id, aug_alts, tm_alts, bins, binned_transcripts):
    """
    bins analyzed candidates based on tier
    """
    for b, vals in itertools.izip(bins, [[winner_id], aug_alts, tm_alts]):
        for v in vals:
            binned_transcripts[b].append(v)


def bin_transcripts(reverse_name_map, stats_dict, ok_ids, ens_ids):
    binned_transcripts = {"bestOk": [], "augAltOk": [], "tmAltOk": [], "bestNotOk": [], "augAltNotOk": [], 
                          "tmAltNotOk": [], "fail": [], "discarded": [], "tieIds": set()}
    for ens_id in ens_ids:
        aln_ids = reverse_name_map[ens_id]
        ok_candidates, not_ok_candidates, discarded = find_ok_not_ok_candidates(stats_dict, aln_ids, ok_ids,
                                                                                discard_cov_cutoff=0.50,
                                                                                filter_cov_cutoff=0.80)
        binned_transcripts["discarded"].extend(discarded)
        if len(ok_candidates) == len(not_ok_candidates) == 0:
            binned_transcripts["fail"].append(ens_id)
        elif len(ok_candidates) > 0:
            ok_id, aug_alt, tm_alt, is_tie = analyze_candidates(ok_candidates)
            bin_candidates(ok_id, aug_alt, tm_alt, ["bestOk", "augAltOk", "tmAltOk"], binned_transcripts)
            if is_tie:
                binned_transcripts["tieIds"].add(strip_alignment_numbers(ok_id))
        else:
            not_ok_id, aug_alt, tm_alt, is_tie = analyze_candidates(not_ok_candidates)
            bin_candidates(not_ok_id, aug_alt, tm_alt, ["bestNotOk", "augAltNotOk", "tmAltNotOk"], binned_transcripts)
            if is_tie:
                binned_transcripts["tieIds"].add(strip_alignment_numbers(not_ok_id))
    return binned_transcripts


def fix_gene_pred(gp, gene_map):
    """
    These genePreds have a few problems. First, the alignment numbers must be removed. Second, we want to fix
    the name2 field to be the gene name. Third, we want to set the unique ID field. Finally, we want to sort the whole 
    thing by genomic coordinates.
    """
    gp = sorted([x.split("\t") for x in gp], key=lambda x: [x[1], x[3]])
    fixed = []
    for x in gp:
        x[10] = x[0]  # use unique Aug/TM ID as unique identifier
        tx_id = strip_alignment_numbers(x[0])
        x[0] = tx_id
        gene_id = gene_map[tx_id]
        x[11] = gene_id
        fixed.append(x)
    return ["\t".join(x) for x in fixed]


def write_gps(binned_transcripts, gps, out_dir, genome, biotype, gene_map):
    p = os.path.join(out_dir, genome, biotype)
    mkdir_p(p)
    for b in ["bestOk", "augAltOk", "tmAltOk", "bestNotOk", "augAltNotOk", "tmAltNotOk"]:
        gp = [gps[aln_id] for aln_id in binned_transcripts[b]]
        fixed_gp = fix_gene_pred(gp, gene_map)
        with open(os.path.join(p, genome + "_" + b + ".gp"), "w") as outf:
            for x in fixed_gp:
                outf.write(x)
    for b in ["fail", "discarded"]:
        with open(os.path.join(p, genome + "_" + b + ".txt"), "w") as outf:
            for aln_id in binned_transcripts[b]:
                outf.write(aln_id + "\n")


def consensus_gene_set(binned_transcripts, stats_dict, gps, gene_map, gene_biotype_map, gene_cov_cutoff=0.20):
    """
    Builds the consensus gene set. For each transcript that has a best OK/not OK (passes coverage filter),
    report it. For the remaining transcripts, determine the set of genes they come from. Find anything above
    gene_cov_cutoff to be the one transcript to best represent this gene.
    """
    consensus = []
    best_ids = set()
    for b in ["bestOk", "bestNotOk"]:
        best_ids |= set(binned_transcripts[b])
        for aln_id in binned_transcripts[b]:
            consensus.append(gps[aln_id])
    discarded_genes = defaultdict(list)
    genes_we_have = {gene_map[strip_alignment_numbers(x)] for x in best_ids}
    for x in binned_transcripts["discarded"]:
        gene = gene_map[strip_alignment_numbers(x)]
        if gene not in genes_we_have:
            discarded_genes[gene].append(x)
    for gene, transcripts in discarded_genes.iteritems():
        stats = [stats_dict[x] for x in transcripts if stats_dict[x][1] != None and stats_dict[x][2] > gene_cov_cutoff]
        if len(stats) > 0:
            best = find_best_aln(stats)[0]
            consensus.append(gps[best])
    return consensus


def write_consensus(consensus, gene_map, consensus_path):
    consensus = fix_gene_pred(consensus, gene_map)
    with open(consensus_path, "w") as outf:
        for x in consensus:
            outf.write(x)


def make_tx_counts_dict(binned_transcripts, filter_set=set()):
    """
    Makes a counts dictionary from binned_transcripts.
    """
    counts = OrderedDict()
    for key in ['sameOk', 'augOk', 'tmOk', 'sameNotOk', 'augNotOk', 'tmNotOk',  'fail']:
        counts[key] = 0
    tie_ids = binned_transcripts["tieIds"]
    for x in binned_transcripts["bestOk"]:
        aln_id = strip_alignment_numbers(x)
        if aln_id not in filter_set:
            continue
        elif aln_id in tie_ids:
            counts["sameOk"] += 1
        elif 'aug' in x:
            counts["augOk"] += 1
        else:
            counts["tmOk"] += 1
    for x in binned_transcripts["bestNotOk"]:
        aln_id = strip_alignment_numbers(x)
        if aln_id not in filter_set:
            continue
        elif aln_id in tie_ids:
            counts["sameNotOk"] += 1
        elif 'aug' in x:
            counts["augNotOk"] += 1
        else:
            counts["tmNotOk"] += 1
    counts["fail"] = len(binned_transcripts["fail"])
    return counts


def main():
    args = parse_args()
    con, cur = sql_lib.attach_databases(args.compAnnPath, has_augustus=True)
    biotypes = get_all_biotypes(args.attributePath)
    gene_map = get_gene_map(args.attributePath)
    gene_biotype_map = get_gene_biotype_map(args.attributePath)
    chr_y_ids = gp_chrom_filter(args.compGp)
    consensus_base_path = os.path.join(args.outDir, "geneSets")
    mkdir_p(consensus_base_path)
    raw_bins_base_path = os.path.join(args.outDir, "binnedTranscripts")
    mkdir_p(raw_bins_base_path)
    binned_transcript_holder = defaultdict(dict)  # save all bins to make some plots at the end
    consensus = []
    coding_ok = get_all_ok(cur, args.genome)
    noncoding_ok = transmap_ok(cur, args.genome, coding=False)
    gps = load_gps([args.tmGp, args.augGp])
    for biotype in biotypes:
        ens_ids = get_all_ids(args.attributePath, biotype=biotype) - chr_y_ids  # filter out chrY
        reverse_name_map = get_reverse_name_map(cur, args.genome, whitelist=ens_ids, has_augustus=True)
        stats_dict = merge_stats(cur, args.genome)
        if biotype == "protein_coding": 
            binned_transcripts = bin_transcripts(reverse_name_map, stats_dict, coding_ok, ens_ids)
        else:
            binned_transcripts = bin_transcripts(reverse_name_map, stats_dict, noncoding_ok, ens_ids)
        consensus.extend(consensus_gene_set(binned_transcripts, stats_dict, gps, gene_map, gene_biotype_map))
        write_gps(binned_transcripts, gps, raw_bins_base_path, args.genome, biotype, gene_map)
        binned_transcript_holder[biotype] = binned_transcripts
    consensus_path = os.path.join(consensus_base_path, args.genome + "_consensusGeneSet.gp")
    write_consensus(consensus, gene_map, consensus_path)
    with open(args.binnedTranscriptPath, 'w') as outf:
        pickle.dump(binned_transcript_holder, outf) 


if __name__ == "__main__":
    main()