import os
import itertools
import argparse
import sqlite3 as sql
from collections import defaultdict
from lib.psl_lib import removeAlignmentNumber, removeAugustusAlignmentNumber
from lib.sqlite_lib import attachDatabase
from lib.general_lib import mkdir_p
from sonLib.bioio import system


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--compAnnPath", required=True)
    parser.add_argument("--statsDir", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--attributePath", required=True)
    parser.add_argument("--augGps", nargs="+", required=True)
    parser.add_argument("--tmGps", nargs="+", required=True)
    return parser.parse_args()


# these classifiers define OK for coding transcripts
tm_coding_classifiers = ["CodingInsertions", "CodingDeletions", "StartOutOfFrame", "FrameShift", 
                         "AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "BadFrame", "BeginStart",
                         "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice", 
                         "EndStop", "InFrameStop", "ShortCds", "UnknownBases", "AlignmentAbutsUnknownBases"]

# these classifiers define OK for non-coding transcripts
tm_noncoding_classifiers = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UtrUnknownSplice",
                            "UtrGap", "UnknownGap", "UnknownBases", "AlignmentAbutsUnknownBases"]


def skip_header(path):
    f_h = open(path)
    _ = f_h.next()
    return f_h


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


def get_all_ids(attr_path, biotype=None):
    """
    returns the set of ensembl IDs in the entire Gencode database pulled from the attribute
    """
    if biotype is None:
        return {x.split()[3] for x in skip_header(attr_path)}
    else:
        return {x.split()[3] for x in skip_header(attr_path) if x.split()[4] == biotype}


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
        ens_id = removeAlignmentNumber(removeAugustusAlignmentNumber(aln_id))
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


def find_t1_t3_candidates(stats_dict, ids, ok_ids, discard_cov_cutoff, filter_cov_cutoff):
    """
    Returns a list of candidate transcripts for tier1 and tier3. All transcripts must have discard_cov_cutoff
    coverage or they are discarded. In order to be a t1 candidate, a transcript must be 'OK' and have coverage above
    filter_cov_cutoff.
    """
    t1_candidates = []
    t3_candidates = []
    discarded = []
    for aln_id in ids:
        if stats_dict[aln_id][2] < discard_cov_cutoff:
            discarded.append(stats_dict[aln_id][0])
        elif aln_id in ok_ids and stats_dict[aln_id][2] >= filter_cov_cutoff:
            t1_candidates.append(stats_dict[aln_id])
        else:
            t3_candidates.append(stats_dict[aln_id])
    return t1_candidates, t3_candidates, discarded


def find_best_aln(stats):
    """
    Takes the list of t1/t3 candidates from find_t1_t3_candidates and returns the best alignment
    """
    return sorted(stats, key=lambda x: -x[1])[0][0]


def split_alternatives(stats, best_id):
    """
    Takes the list of t1/t3 candidates and splits them up by augustus/transMap, filtering out the winner
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
    or transMap
    """
    winner_id = find_best_aln(candidates)
    aug_alts, tm_alts = split_alternatives(candidates, winner_id)
    return winner_id, aug_alts, tm_alts


def bin_candidates(winner_id, aug_alts, tm_alts, bins, binned_transcripts):
    """
    bins analyzed candidates based on tier
    """
    for b, vals in itertools.izip(bins, [[winner_id], aug_alts, tm_alts]):
        for v in vals:
            binned_transcripts[b].append(v)


def bin_transcripts(reverse_name_map, stats_dict, ok_ids, ens_ids):
    binned_transcripts = {"bestOk": [], "augAltOk": [], "tmAltOk": [], "bestNotOk": [], "augAltNotOk": [], 
                          "tmAltNotOk": [], "fail": [], "discarded": []}
    for ens_id in ens_ids:
        aln_ids = reverse_name_map[ens_id]
        t1_candidates, t3_candidates, discarded = find_t1_t3_candidates(stats_dict, aln_ids, ok_ids,
                                                                        discard_cov_cutoff=0.50, filter_cov_cutoff=0.80)
        binned_transcripts["discarded"].extend(discarded)
        if len(t1_candidates) == len(t3_candidates) == 0:
            binned_transcripts["fail"].append(ens_id)
        elif len(t1_candidates) > 0:
            t1_id, aug_t2, tm_t2 = analyze_candidates(t1_candidates)
            bin_candidates(t1_id, aug_t2, tm_t2, ["bestOk", "augAltOk", "tmAltOk"], binned_transcripts)
        else:
            t3_id, aug_t4, tm_t4 = analyze_candidates(t3_candidates)
            bin_candidates(t3_id, aug_t4, tm_t4, ["bestNotOk", "augAltNotOk", "tmAltNotOk"], binned_transcripts)
    return binned_transcripts


def load_gps(gp_paths):
    return {l.split()[0]: l for p in gp_paths for l in open(p)}


def fix_gene_pred(gp, gene_map):
    """
    These genePreds have a few problems. First, the alignment numbers must be removed. Second, we want to fix
    the name2 field to be the gene name. Finally, we want to sort the whole thing by genomic coordinates.
    """
    gp = sorted([x.split("\t") for x in gp], key=lambda x: [x[1], x[3]])
    fixed = []
    for i, x in enumerate(gp):
        tx_id = removeAlignmentNumber(removeAugustusAlignmentNumber(x[0]))
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
    return {gene_map[removeAlignmentNumber(removeAugustusAlignmentNumber(x))] for x in not_ok_gps}


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
        discarded_genes[gene_map[removeAlignmentNumber(removeAugustusAlignmentNumber(x))]].append(x)
    for gene, transcripts in discarded_genes.iteritems():
        stats = [stats_dict[x] for x in transcripts if stats_dict[x][1] != None and stats_dict[x][2] > gene_cov_cutoff]
        if len(stats) > 0:
            best = find_best_aln(stats)
            consensus.append(gps[best])
    return consensus


def write_consensus(consensus, gene_map, consensus_path):
    consensus = fix_gene_pred(consensus, gene_map)
    with open(consensus_path, "w") as outf:
        for x in consensus:
            outf.write(x)


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
    raw_bins_base_path = os.path.join(args.outDir, "binnedTranscripts")
    mkdir_p(raw_bins_base_path)
    for genome, tm_gp, aug_gp in itertools.izip(sorted_genomes, sorted_tm_gps, sorted_aug_gps):
        consensus = []
        assert genome in tm_gp and genome in aug_gp
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
        consensus_path = os.path.join(consensus_base_path, genome + "consensusGeneSet.gp")
        write_consensus(consensus, gene_map, consensus_path)

if __name__ == "__main__":
    main()