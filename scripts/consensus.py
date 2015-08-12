import os
import itertools
import argparse
import sqlite3 as sql
from lib.psl_lib import removeAlignmentNumber, removeAugustusAlignmentNumber
from lib.sqlite_lib import attachDatabase
from lib.general_lib import mkdir_p


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


def transmap_ok(cur, genome):
    """
    Finds all aIds which are 'OK' based on the classifyFields below
    """
    classifyFields = ["CodingInsertions", "CodingDeletions", "CodingDeletions", "StartOutOfFrame", "FrameShift", 
                      "AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "BadFrame", "BeginStart",
                      "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice", 
                      "EndStop", "InFrameStop", "ShortCds", "UnknownBases", "AlignmentAbutsUnknownBases"]
    cmd = """SELECT main.'{0}'.'AlignmentId' FROM main.'{0}' WHERE (""".format(genome)
    for col in classifyFields[:-1]:
        cmd += " main.'{}'.'{}' = ? {}".format(genome, col, "AND")
    cmd += " main.'{}'.'{}' = ?)".format(genome, classifyFields[-1])
    vals = [0] * len(classifyFields)
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


def get_all_ok(cur, genome):
    """
    Adapter function to return the combined sets of ok from augustus and transmap
    """
    return augustus_ok(cur, genome) | transmap_ok(cur, genome)


def get_all_ids(attr_path, biotype=None):
    """
    returns the set of ensembl IDs in the entire Gencode database pulled from the attribute
    """
    if biotype is None:
        return {x.split()[3] for x in open(attr_path)}
    else:
        return {x.split()[3] for x in open(attr_path) if x.split()[4] == biotype}


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
    binned_transcripts = {"T1": [], "T2A": [], "T2T": [], "T3": [], "T4A": [], "T4T": [], "fail": [], "discarded": []}
    for ens_id in ens_ids:
        aln_ids = reverse_name_map[ens_id]
        t1_candidates, t3_candidates, discarded = find_t1_t3_candidates(stats_dict, aln_ids, ok_ids,
                                                                        discard_cov_cutoff=0.50, filter_cov_cutoff=0.80)
        binned_transcripts["discarded"].extend(discarded)
        if len(t1_candidates) == len(t3_candidates) == 0:
            binned_transcripts["fail"].append(ens_id)
        elif len(t1_candidates) > 0:
            t1_id, aug_t2, tm_t2 = analyze_candidates(t1_candidates)
            bin_candidates(t1_id, aug_t2, tm_t2, ["T1", "T2A", "T2T"], binned_transcripts)
        else:
            t3_id, aug_t4, tm_t4 = analyze_candidates(t3_candidates)
            bin_candidates(t3_id, aug_t4, tm_t4, ["T3", "T4A", "T4T"], binned_transcripts)
    return binned_transcripts


def load_gps(gp_paths):
    return {l.split()[0]: l for p in gp_paths for l in open(p)}


def write_gps(binned_transcripts, gps, out_dir, genome, filter_set=None):
    p = os.path.join(out_dir, genome)
    mkdir_p(p)
    for b in ["T1", "T2A", "T2T", "T3", "T4A", "T4T"]:
        with open(os.path.join(p, genome + "_Tier" + b + ".gp"), "w") as outf:
            for aln_id in binned_transcripts[b]:
                gp = gps[aln_id]
                if filter_set is None:
                    outf.write(gp)
                elif removeAlignmentNumber(removeAugustusAlignmentNumber(aln_id)) in filter_set:
                    outf.write(gp)
    for b in ["fail", "discarded"]:
        with open(os.path.join(p, genome + "_" + b + ".txt"), "w") as outf:
            for aln_id in binned_transcripts[b]:
                if filter_set is None:
                    outf.write(aln_id + "\n")
                elif removeAlignmentNumber(removeAugustusAlignmentNumber(aln_id)) in filter_set:
                    outf.write(aln_id + "\n")


def driver_function(genome, tm_gp, aug_gp, con, cur, ens_ids, coding_ids, stats_dir, out_dir):
    reverse_name_map = get_reverse_name_map(cur, genome, ens_ids)
    ok_ids = get_all_ok(cur, genome)
    stats_dict = merge_stats(cur, stats_dir, genome)
    binned_transcripts = bin_transcripts(reverse_name_map, stats_dict, ok_ids, ens_ids)
    gps = load_gps([tm_gp, aug_gp])
    all_path = os.path.join(out_dir, "all_transcripts")
    write_gps(binned_transcripts, gps, all_path, genome)
    coding_path = os.path.join(out_dir, "coding_transcripts")
    write_gps(binned_transcripts, gps, coding_path, genome, filter_set=coding_ids) 


def main():
    args = parse_args()
    con, cur = attach_databases(args.compAnnPath)
    ens_ids = get_all_ids(args.attributePath)
    coding_ids = get_all_ids(args.attributePath, biotype="protein_coding")
    sorted_genomes = sorted(args.genomes)
    sorted_tm_gps = sorted(args.tmGps)
    sorted_aug_gps = sorted(args.augGps)
    for genome, tm_gp, aug_gp in itertools.izip(sorted_genomes, sorted_tm_gps, sorted_aug_gps):
        assert genome in tm_gp and genome in aug_gp
        driver_function(genome, tm_gp, aug_gp, con, cur, ens_ids, coding_ids, args.statsDir,
                        args.outDir)


if __name__ == "__main__":
    main()