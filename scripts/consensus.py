import os
import itertools
from lib.psl_lib import *  # need to convert names
from lib.sqlite_lib import *
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.lines as lines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.backends.backend_pdf as pltBack
import numpy as np
from collections import OrderedDict
from scripts.coverage_identity_ok_plots import *
from scripts.alignAugustus import mkdir_p

def transmap_ok(cur, genome):
    """
    Finds all aIds which are 'OK' based on the classifyFields below
    """
    classifyFields = ["CodingInsertions","CodingDeletions", "CodingDeletions", "StartOutOfFrame", "FrameShift", "AlignmentAbutsLeft", "AlignmentAbutsRight", 
                      "AlignmentPartialMap", "BadFrame", "BeginStart", "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice", 
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
    classifyFields = ["AugustusSameStartStop", "AugustusExonGain", "AugustusExonLoss", "AugustusParalogy",
                      "AugustusNotSimilarExonBoundaries"]
    cmd = """SELECT augustus.'{0}'.'AlignmentId' FROM augustus.'{0}' WHERE (""".format(genome)
    for col in classifyFields[:-1]:
        cmd += " augustus.'{}'.'{}' = ? {}".format(genome, col, "AND")
    cmd += " augustus.'{}'.'{}' = ?)".format(genome, classifyFields[-1])
    vals = [0] * len(classifyFields)
    return {x[0] for x in cur.execute(cmd, vals).fetchall()}


def get_all_ids(attr_path):
    """
    returns the set of ensembl IDs in the entire Gencode database pulled from the attribute
    """
    return {x.split()[0] for x in open(attr_path)}


def reverse_name_map(cur, genome, ids):
    """
    creates a dictionary mapping each Gencode ID to all IDs produced by Augustus and transMap
    """
    reverse_map = {x: [[], []] for x in ids}
    base_cmd = "SELECT {0}.'{1}'.'AlignmentId' FROM {0}.'{1}'"
    aug_cmd = base_cmd.format("augustus", genome)
    tm_cmd = base_cmd.format("main", genome)
    aug_r = cur.execute(aug_cmd).fetchall()
    tm_r = cur.execute(tm_cmd).fetchall()
    for aln_id in aug_r:
        aln_id = aln_id[0]
        ens_id = removeAlignmentNumber(removeAugustusAlignmentNumber(aln_id))
        reverse_map[ens_id][0].append(aln_id)
    for aln_id in tm_r:
        aln_id = aln_id[0]
        ens_id = removeAlignmentNumber(removeAugustusAlignmentNumber(aln_id))
        reverse_map[ens_id][1].append(aln_id)
    return reverse_map


def get_tm_stats(cur, genome):
    """
    Pulls the alignment metrics from the attributes database
    """
    cmd = "SELECT AlignmentId, AlignmentIdentity, AlignmentCoverage FROM attributes.'{}' WHERE AlignmentIdentity " \
          "IS NOT NULL AND AlignmentCoverage IS NOT NULL".format(genome)
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


def attach_databases(base_path):
    """
    Attaches all of the databases. Expects base_path to be the path that comparativeAnnotator wrote to.
    """
    con = sql.connect(os.path.join(base_path, "classify.db"))
    cur = con.cursor()
    attachDatabase(con, os.path.join(base_path, "augustusClassify.db"), "augustus")
    attachDatabase(con, os.path.join(base_path, "attributes.db"), "attributes")
    return con, cur


def filter_transcripts(t, cov_cutoff):
    """
    Filters a list of transcripts in the form [id, ident, cov] based on cov_cutoff
    and returns this list sorted by identity
    """
    return sorted([x for x in t if x[2] >= cov_cutoff], key=lambda x: -x[1])


def find_t1_t3(aug_stats, tm_stats, cov_cutoff):
    """
    Uses filter_transcripts() to find the tier1/tier3 transcript for this ens_id
    NOTE: due to the sorting, this will always pick the augustus transcript in a tie.
    This may inflate how much augustus is supposedly improving things.
    """
    t = filter_transcripts(itertools.chain(aug_stats, tm_stats), cov_cutoff)
    if len(t) > 0:
        return t[0][0]
    else:
        return None


def find_t2_t4(stats, best_transcript, cov_cutoff):
    """
    Uses filter_transcripts() to find the tier2/tier4 transcripts for this ens_id
    Transcripts which fail the coverage cutoff will be reported separately to put in tier4 later
    """
    t = filter_transcripts(stats, cov_cutoff)
    return [x[0] for x in t if x[0] != best_transcript], [x[0] for x in stats if x not in t]


def bin_transcripts(reverse_map, aug_ok_ids, tm_ok_ids, aug_stats, tm_stats, ids):
    binned_transcripts = {"T1": [], "T2A": [], "T2T": [], "T3": [], "T4A": [], "T4T": [], "fail": []}
    for ens_id in ids:
        aug_ids, tm_ids = reverse_map[ens_id]
        if len(aug_ids) == len(tm_ids) == 0:
            binned_transcripts["fail"].append(ens_id)
            continue  # this is a tier5 ens_id - we have nothing for it
        aug_ok_stats = [aug_stats[x] for x in aug_ids if x in aug_ok_ids]
        tm_ok_stats = [tm_stats[x] for x in tm_ids if x in tm_ok_ids]
        t1_transcript = find_t1_t3(aug_ok_stats, tm_ok_stats, cov_cutoff=0.85)
        if t1_transcript is not None:
            binned_transcripts["T1"].append(t1_transcript)
        t2a_transcripts, t2a_fail = find_t2_t4(aug_ok_stats, t1_transcript, cov_cutoff=0.50)
        t2t_transcripts, t2t_fail = find_t2_t4(tm_ok_stats, t1_transcript, cov_cutoff=0.50)
        binned_transcripts["T2A"].extend(t2a_transcripts)
        binned_transcripts["T2T"].extend(t2t_transcripts)
        aug_not_ok_stats = [aug_stats[x] for x in aug_ids if x not in aug_ok_ids or x in t2a_fail]
        tm_not_ok_stats = [tm_stats[x] for x in tm_ids if x not in tm_ok_ids or x in t2t_fail]
        if t1_transcript is None:
            t3_transcript = find_t1_t3(aug_not_ok_stats, tm_not_ok_stats, cov_cutoff=0.50)
            if t3_transcript is not None:
                binned_transcripts["T3"].append(t3_transcript)
            else:
                binned_transcripts["fail"].append(ens_id)
        else:
            t3_transcript = None
        t4a_transcripts, t4a_fail = find_t2_t4(aug_not_ok_stats, t3_transcript, cov_cutoff=0.50)
        t4t_transcripts, t4t_fail = find_t2_t4(tm_not_ok_stats, t3_transcript, cov_cutoff=0.50)
        binned_transcripts["T4A"].extend(t4a_transcripts)
        binned_transcripts["T4T"].extend(t4t_transcripts)
    return binned_transcripts
# error is here: tm_stats[x] why is x not in tm_stats?: