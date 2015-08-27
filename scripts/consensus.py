import os
import itertools
import argparse
import sqlite3 as sql
from collections import defaultdict, Counter, OrderedDict
from lib.psl_lib import removeAlignmentNumber, removeAugustusAlignmentNumber
from lib.sqlite_lib import attachDatabase
from lib.general_lib import mkdir_p, DefaultOrderedDict
from sonLib.bioio import system

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.lines as lines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.backends.backend_pdf as pltBack
import numpy as np
from scripts.coverage_identity_ok_plots import init_image, establish_axes, plot_bars


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
    return parser.parse_args()


# these classifiers define OK for coding transcripts
tm_coding_classifiers = ["CodingInsertions", "CodingDeletions", "StartOutOfFrame", "FrameShift", 
                         "AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "BadFrame", "BeginStart",
                         "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice", 
                         "EndStop", "InFrameStop", "ShortCds", "UnknownBases", "AlignmentAbutsUnknownBases"]

# these classifiers define OK for non-coding transcripts
tm_noncoding_classifiers = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UtrUnknownSplice",
                            "UtrGap", "UnknownGap", "UnknownBases", "AlignmentAbutsUnknownBases"]

# used for the plots
width=8.0
height=4.0
bar_width=0.4

def skip_header(path):
    """
    The attributes file produced by the pipeline has a header. Skip it. Return a open file handle pointing to line 2.
    """
    f_h = open(path)
    _ = f_h.next()
    return f_h


def strip_alignment_numbers(aln_id):
    """
    Convenience function for stripping both Augustus and transMap alignment IDs from a aln_id
    """
    return removeAlignmentNumber(removeAugustusAlignmentNumber(aln_id))


def get_all_biotypes(attr_path):
    """
    Returns all biotypes in the attribute database.
    """
    return {x.split()[4] for x in skip_header(attr_path)}


def transmap_ok(cur, genome, classify_fields):
    """
    Finds all aIds which are 'OK' based on the classifyFields below
    """
    cmd = """SELECT main.'{0}'.'AlignmentId' FROM main.'{0}' WHERE (""".format(genome)
    for col in classify_fields[:-1]:
        cmd += " main.'{}'.'{}' = ? {}".format(genome, col, "AND")
    cmd += " main.'{}'.'{}' = ?)".format(genome, classify_fields[-1])
    vals = [0] * len(classify_fields)
    return {x[0] for x in cur.execute(cmd, vals).fetchall()}


def augustus_ok(cur, genome):
    """
    Finds all aug_aIds which are 'OK' as defined by the fields in classifyFields
    """
    classifyFields = ['AugustusParalogy', 'AugustusExonGain', 'AugustusExonLoss', 'AugustusNotSameStrand', 
                      'AugustusNotSameStartStop', 'AugustusNotSimilarTerminalExonBoundaries', 
                      'AugustusNotSimilarInternalExonBoundaries']
    cmd = """SELECT augustus.'{0}'.'AlignmentId' FROM augustus.'{0}' WHERE (""".format(genome)
    for col in classifyFields[:-1]:
        cmd += " augustus.'{}'.'{}' = ? {}".format(genome, col, "AND")
    cmd += " augustus.'{}'.'{}' = ?)".format(genome, classifyFields[-1])
    vals = [0] * len(classifyFields)
    return {x[0] for x in cur.execute(cmd, vals).fetchall()}


def get_all_ok(cur, genome, tm_classifiers):
    """
    Adapter function to return the combined sets of ok from augustus and transmap
    """
    return augustus_ok(cur, genome) | transmap_ok(cur, genome, tm_classifiers)


def get_all_ids(attr_path, biotype=None, filter_set=set()):
    """
    returns the set of ensembl IDs in the entire Gencode database pulled from the attribute
    """
    if biotype is None:
        return {x.split()[3] for x in skip_header(attr_path) if x not in filter_set}
    else:
        return {x.split()[3] for x in skip_header(attr_path) if x.split()[4] == biotype if x not in filter_set}


def get_gp_ids(gp):
    return {x.split()[0] for x in open(gp)}


def get_reverse_name_map(cur, genome, ids):
    """
    creates a dictionary mapping each Gencode ID to all IDs produced by Augustus and transMap
    """
    reverse_name_map = {x: [] for x in ids}
    base_cmd = "SELECT {0}.'{1}'.'AlignmentId' FROM {0}.'{1}'"
    aug_cmd = base_cmd.format("augustus", genome)
    tm_cmd = base_cmd.format("main", genome)
    aug_r = cur.execute(aug_cmd).fetchall()
    tm_r = cur.execute(tm_cmd).fetchall()
    for aln_id in itertools.chain(aug_r, tm_r):
        aln_id = aln_id[0]
        ens_id = strip_alignment_numbers(aln_id)
        if ens_id in ids:
            reverse_name_map[ens_id].append(aln_id)
    return reverse_name_map


def get_tm_stats(cur, genome):
    """
    Pulls the alignment metrics from the attributes database
    """
    cmd = "SELECT AlignmentId, AlignmentIdentity, AlignmentCoverage FROM attributes.'{}'".format(genome)
    result = cur.execute(cmd).fetchall()
    return {x[0]: x for x in result}


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


def attach_databases(comp_ann_path):
    """
    Attaches all of the databases. Expects comp_ann_path to be the path that comparativeAnnotator wrote to.
    """
    con = sql.connect(os.path.join(comp_ann_path, "classify.db"))
    cur = con.cursor()
    attachDatabase(con, os.path.join(comp_ann_path, "augustusClassify.db"), "augustus")
    attachDatabase(con, os.path.join(comp_ann_path, "attributes.db"), "attributes")
    return con, cur


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


def load_gps(gp_paths):
    return {l.split()[0]: l for p in gp_paths for l in open(p)}


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


def find_not_ok_genes(gene_map, not_ok_gps):
    """
    Returns the set of (not failed/discarded) genes for which we have no OK transcripts.
    """
    return {gene_map[strip_alignment_numbers(x)] for x in not_ok_gps}


def get_gene_biotype_map(attr_path):
    """
    Returns a dictionary mapping all gene IDs to their respective biotypes
    """
    return {x.split()[0]: x.split()[2] for x in skip_header(attr_path)}


def get_gene_map(attr_path):
    """
    Returns a dictionary mapping all transcript IDs to their respective gene IDs
    """
    return {x.split()[3]: x.split()[0] for x in skip_header(attr_path)}


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


def make_counts_frequency(counts):
    """
    Convenience function that takes a dict and turns the values into a proportion of the total.
    Returns a list of lists [[name, percent]]
    """
    tot = sum(counts.values())
    for key, val in counts.iteritems():
        counts[key] = 1.0 * val / tot
    return list(counts.iteritems())


def make_counts_dict(binned_transcripts, filter_set=set()):
    """
    Makes a counts dictionary from binned_transcripts.
    """
    counts = OrderedDict()
    for key in ['augOk', 'tmOk', 'sameOk', 'augNotOk', 'tmNotOk', 'sameNotOk', 'fail']:
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
    Defines a fixed order of genomes based on the most OK protein coding
    """
    ok_counts = []
    for genome, binned_transcripts in binned_transcript_holder.iteritems():
        num_ok = len(binned_transcripts["protein_coding"]["bestOk"])
        ok_counts.append([genome, 1.0 * num_ok / len(ens_ids)])
    genome_order = sorted(ok_counts, key=lambda x: -x[1])
    return zip(*genome_order)[0]


def barplot(results, color_palette, out_path, file_name, title_string, categories, border=True, has_legend=True):
    """
    Boilerplate code that will produce a barplot.
    """
    fig, pdf = init_image(out_path, file_name, width, height)
    ax = establish_axes(fig, width, height, border, has_legend)
    plt.text(0.5, 1.08, title_string, horizontalalignment='center', fontsize=12, transform=ax.transAxes)
    ax.set_ylabel("Proportion of transcripts")
    ax.set_ylim([0, 1.0])
    plt.tick_params(axis='y', labelsize=8)
    plt.tick_params(axis='x', labelsize=9)
    ax.yaxis.set_ticks(np.arange(0.0, 101.0, 10.0) / 100.0)
    ax.yaxis.set_ticklabels([str(x) + "%" for x in range(0, 101, 10)])
    ax.xaxis.set_ticks(np.arange(0, len(results)) + bar_width / 2.0)
    ax.xaxis.set_ticklabels(zip(*results)[0], rotation=55)
    bars = plot_bars(ax, zip(*results)[1], bar_width, color_palette=color_palette)
    legend = fig.legend([x[0] for x in bars[::-1]], categories[::-1], bbox_to_anchor=(1,0.8), fontsize=11, frameon=True, 
                        title="Category")
    fig.savefig(pdf, format='pdf')
    pdf.close()


def make_coding_transcript_plot(binned_transcript_holder, out_path, out_name, genome_order, ens_ids, title_string):
    protein_coding_palette = ["#df65b0", "#dd1c77", "#980043", "#a1dab4", "#41b6c4", "#2c7fb8", "#252525"]
    coding_metrics = OrderedDict()
    for g in genome_order:
        bins = binned_transcript_holder[g]['protein_coding']
        coding_metrics[g] = make_counts_frequency(make_counts_dict(bins, filter_set=ens_ids))
    categories = zip(*coding_metrics[g])[0]
    results = [[g, zip(*coding_metrics[g])[1]] for g in coding_metrics]
    barplot(results, protein_coding_palette, out_path, out_name, title_string, categories)


base_title_string = "Proportion of {:,} {}\nOK / notOK In Target Genomes"
title_string_dict = {"Comp": "Protein-Coding Transcripts in GencodeCompVM4", 
                     "Basic": "Protein-Coding Transcripts in GencodeBasicVM4", 
                     "Complement": "Protein-coding Transcripts\nin GencodeCompVM4 and NOT in GencodeBasicVM4"}
file_name_dict = {"Comp": "protein_coding_comprehensive", "Basic": "protein_coding_basic", 
                  "Complement": "protein_coding_complement"}


def make_coding_plots(binned_transcript_holder, out_path, comp_gp, basic_gp, attr_path):
    comp_ids = get_gp_ids(comp_gp)
    basic_ids = get_gp_ids(basic_gp)
    coding_ids = get_all_ids(attr_path, biotype="protein_coding")
    basic_coding = basic_ids & coding_ids
    comp_coding = comp_ids & coding_ids
    complement_coding = comp_coding - basic_coding
    genome_order = find_genome_order(binned_transcript_holder, comp_coding)
    for cat, ids in zip(*[["Comp", "Basic", "Complement"], [comp_coding, basic_coding, complement_coding]]):
        title_string = base_title_string.format(len(ids), title_string_dict[cat])
        out_name = file_name_dict[cat]
        make_coding_transcript_plot(binned_transcript_holder, out_path, out_name, genome_order, ids, cat)


def main():
    args = parse_args()
    con, cur = attach_databases(args.compAnnPath)
    biotypes = get_all_biotypes(args.attributePath)
    gene_map = get_gene_map(args.attributePath)
    gene_biotype_map = get_gene_biotype_map(args.attributePath)
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
            ens_ids = get_all_ids(args.attributePath, biotype=biotype)
            reverse_name_map = get_reverse_name_map(cur, genome, ens_ids)
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
    make_coding_plots(binned_transcript_holder, plots_path, args.compGp, args.basicGp, args.attributePath)

if __name__ == "__main__":
    main()