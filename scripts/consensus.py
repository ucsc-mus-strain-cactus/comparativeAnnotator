import argparse
from scripts.plot_functions import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--statsDir", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--attributePath", required=True)
    parser.add_argument("--augGps", nargs="+", required=True)
    parser.add_argument("--tmGps", nargs="+", required=True)
    parser.add_argument("--compGp", required=True)
    parser.add_argument("--basicGp", required=True)
    parser.add_argument("--plotBiotypes", default=["protein_coding", "lincRNA", "miRNA", "snoRNA"])
    return parser.parse_args()


def get_aug_stats(stats_dir, genome):
    """
    Pulls the alignment metrics from the output of alignAugustus.py
    """
    aln_stats = {}
    with open(os.path.join(stats_dir, genome + ".stats")) as f:
        for l in f:
            aug_aId, ident, cov = l.split()
            aln_stats[aug_aId] = [aug_aId] + map(float, [ident, cov])
    return aln_stats


def merge_stats(cur, stats_dir, genome):
    """
    Adapter function to return the combination of the stats dicts produced for augustus and transMap
    """
    tm_stats = get_tm_stats(cur, genome)
    aug_stats = get_aug_stats(stats_dir, genome)
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


def analyze_candidates(candidates):
    """
    Analyzes candidate transcripts, finding the winner and splitting the alternatives based on coming from augustus
    or transMap. If there are multiple equal winners from both transMap and Augustus, reports so.
    """
    winner_ids = find_best_aln(candidates)
    if len(winner_ids) > 1:
        is_tie = True
    else:
        is_tie = False
    winner_id = winner_ids[0]
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
    for i, x in enumerate(gp):
        tx_id = strip_alignment_numbers(x[0])
        x[0] = tx_id
        gene_id = gene_map[tx_id]
        x[10] = str(i)
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
    for b in ["bestOk", "bestNotOk"]:
        for aln_id in binned_transcripts[b]:
            consensus.append(gps[aln_id])
    discarded_genes = defaultdict(list)
    for x in binned_transcripts["discarded"]:
        discarded_genes[gene_map[strip_alignment_numbers(x)]].append(x)
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


def find_genome_order(binned_transcript_holder, ens_ids):
    """
    Defines a fixed order of genomes based on the most OK protein coding.
    This is deprecated in favor of using a hard coded order provided by Joel.
    """
    ok_counts = []
    for genome, binned_transcripts in binned_transcript_holder.iteritems():
        num_ok = len(binned_transcripts["protein_coding"]["bestOk"])
        ok_counts.append([genome, 1.0 * num_ok / len(ens_ids)])
    genome_order = sorted(ok_counts, key=lambda x: -x[1])
    return zip(*genome_order)[0]


def make_coding_transcript_plot(binned_transcript_holder, out_path, out_name, genome_order, ens_ids, title_string):
    coding_metrics = OrderedDict()
    for g in genome_order:
        bins = binned_transcript_holder[g]['protein_coding']
        coding_metrics[g] = make_counts_frequency(make_tx_counts_dict(bins, filter_set=ens_ids))
    categories = zip(*coding_metrics[g])[0]
    results = [[g, zip(*coding_metrics[g])[1]] for g in coding_metrics]
    stacked_barplot(results, categories, out_path, out_name, title_string, color_palette=paired_palette)


base_title_string = "Proportion of {:,} {}\nOK / not OK In Target Genomes"
title_string_dict = {"Comp": "Protein-coding {} in GencodeCompVM4", 
                     "Basic": "Protein-coding {} in GencodeBasicVM4", 
                     "Complement": "Protein-coding {}\nin GencodeCompVM4 and NOT in GencodeBasicVM4"}
file_name_dict = {"Comp": "protein_coding_comprehensive", "Basic": "protein_coding_basic", 
                  "Complement": "protein_coding_complement"}


def make_coding_transcript_plots(binned_transcript_holder, out_path, comp_gp, basic_gp, attr_path):
    comp_ids = get_gp_ids(comp_gp)
    basic_ids = get_gp_ids(basic_gp)
    coding_ids = get_all_ids(attr_path, biotype="protein_coding")
    basic_coding = basic_ids & coding_ids
    comp_coding = comp_ids & coding_ids
    complement_coding = comp_coding - basic_coding
    #genome_order = find_genome_order(binned_transcript_holder, comp_coding)
    genome_order = hard_coded_genome_order
    for cat, ids in zip(*[["Comp", "Basic", "Complement"], [comp_coding, basic_coding, complement_coding]]):
        title_string = base_title_string.format(len(ids), title_string_dict[cat].format("Transcript"))
        out_name = file_name_dict[cat]
        make_coding_transcript_plot(binned_transcript_holder, out_path, out_name, genome_order, ids, title_string)


def calculate_gene_ok_metrics(bins, gene_map, gene_ids):
    ok_genes = {gene_map[strip_alignment_numbers(x)] for x in bins["bestOk"] if strip_alignment_numbers(x) in gene_map}
    not_ok_genes = {gene_map[strip_alignment_numbers(x)] for x in bins["bestNotOk"] if strip_alignment_numbers(x) in 
                    gene_map and gene_map[strip_alignment_numbers(x)] not in ok_genes}
    fail_genes = gene_ids - (ok_genes | not_ok_genes)
    od = OrderedDict([["Has OK Tx", len(ok_genes)], ["No OK Tx", len(not_ok_genes)], ["No Tx", len(fail_genes)]])
    return make_counts_frequency(od)


def ok_gene_by_biotype(binned_transcript_holder, out_path, attr_path, gene_map, genome_order, biotype):
    biotype_ids = get_all_ids(attr_path, biotype=biotype, id_type="Gene")
    title_string = "Proportion of {:,} {} genes with at least one OK transcript".format(len(biotype_ids), biotype)
    file_name = "{}_gene".format(biotype)
    metrics = OrderedDict()
    for g in genome_order:
        bins = binned_transcript_holder[g][biotype]
        metrics[g] = calculate_gene_ok_metrics(bins, gene_map, biotype_ids)
    categories = zip(*metrics[g])[0]
    results = [[g, zip(*metrics[g])[1], zip(*metrics[g])[2]] for g in metrics]
    stacked_barplot(results, categories, out_path, file_name, title_string)


def main():
    args = parse_args()
    con, cur = attach_databases(args.compAnnPath, has_augustus=True)
    biotypes = get_all_biotypes(args.attributePath)
    gene_map = get_gene_map(args.attributePath)
    gene_biotype_map = get_gene_biotype_map(args.attributePath)
    chr_y_ids = gp_chrom_filter(args.compGp)
    sorted_genomes = sorted(args.genomes)
    sorted_tm_gps = sorted(args.tmGps)
    sorted_aug_gps = sorted(args.augGps)
    consensus_base_path = os.path.join(args.outDir, "geneSets")
    mkdir_p(consensus_base_path)
    plots_path = os.path.join(args.outDir, "geneSetMetrics")
    mkdir_p(plots_path)
    raw_bins_base_path = os.path.join(args.outDir, "binnedTranscripts")
    mkdir_p(raw_bins_base_path)
    binned_transcript_holder = defaultdict(dict)  # save all bins to make some plots at the end
    for genome, tm_gp, aug_gp in itertools.izip(sorted_genomes, sorted_tm_gps, sorted_aug_gps):
        consensus = []
        assert genome in tm_gp and genome in aug_gp # sanity check that the right genePreds are being used
        coding_ok = get_all_ok(cur, genome, tm_coding_classifiers)
        noncoding_ok = get_all_ok(cur, genome, tm_noncoding_classifiers)
        gps = load_gps([tm_gp, aug_gp])
        for biotype in biotypes:
            ens_ids = get_all_ids(args.attributePath, biotype=biotype) - chr_y_ids  # filter out chrY
            reverse_name_map = get_reverse_name_map(cur, genome, whitelist=ens_ids, has_augustus=True)
            stats_dict = merge_stats(cur, args.statsDir, genome)
            if biotype == "protein_coding": 
                binned_transcripts = bin_transcripts(reverse_name_map, stats_dict, coding_ok, ens_ids)
            else:
                binned_transcripts = bin_transcripts(reverse_name_map, stats_dict, noncoding_ok, ens_ids)
            consensus.extend(consensus_gene_set(binned_transcripts, stats_dict, gps, gene_map, gene_biotype_map))
            write_gps(binned_transcripts, gps, raw_bins_base_path, genome, biotype, gene_map)
            binned_transcript_holder[genome][biotype] = binned_transcripts
        consensus_path = os.path.join(consensus_base_path, genome + "consensusGeneSet.gp")
        write_consensus(consensus, gene_map, consensus_path)
    make_coding_transcript_plots(binned_transcript_holder, plots_path, args.compGp, args.basicGp, args.attributePath)
    for biotype in args.plotBiotypes:
        ok_gene_by_biotype(binned_transcript_holder, plots_path, attr_path, gene_map, genome_order, biotype)


if __name__ == "__main__":
    main()