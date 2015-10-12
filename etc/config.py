"""
This file configures the plotting and queries for the comparativeAnnotator pipeline.
"""

import src.classifiers
import src.augustus_classifiers
import src.alignment_classifiers
from lib.general_lib import classes_in_module, dict_to_named_tuple, merge_dicts

__author__ = "Ian Fiddes"


# Below are global variables shared by plotting scripts

# this genetic distance is from Joel. There are functions to derive it, but this makes it consistent across plots.
hard_coded_genome_order = ['C57B6NJ', 'NZOHlLtJ', '129S1', 'FVBNJ', 'NODShiLtJ', 'LPJ', 'AJ', 'AKRJ', 'BALBcJ', 'DBA2J',
                            'C3HHeJ', 'CBAJ', 'WSBEiJ', 'CASTEiJ', 'PWKPhJ', 'SPRETEiJ', 'CAROLIEiJ', 'PAHARIEiJ']

# these classifiers define Pass for single-genome analysis
ref_classifiers = ["BadFrame", "BeginStart", "EndStop", "CdsGap", "CdsUnknownSplice", "UtrUnknownSplice",
                   "StartOutOfFrame", "SpliceContainsUnknownBases", "InFrameStop", "ShortCds"]

# these classifiers define Pass/Good for coding transcripts
tm_coding_classifiers = ["BadFrame", "BeginStart", "EndStop", "CdsGap", "CdsUnknownSplice", "UtrUnknownSplice",
                         "StartOutOfFrame", "InFrameStop", "ShortCds", "CodingInsertions", "CodingDeletions",
                         "FrameShift", "HasOriginalStart", "HasOriginalStop", "HasOriginalIntrons"]

# these classifiers define Pass/Good for non-coding transcripts
noncoding_classifiers = ["UtrUnknownSplice"]

# these classifiers define Pass for Augustus transcripts
aug_classifiers = ['AugustusParalogy', 'AugustusExonGain', 'AugustusExonLoss', 'AugustusNotSameStrand',
                      'AugustusNotSameStartStop', 'AugustusNotSimilarTerminalExonBoundaries',
                      'AugustusNotSimilarInternalExonBoundaries']

ref_pass = dict_to_named_tuple({"classifiers": ref_classifiers}, "ref_pass")
aug_pass = dict_to_named_tuple({"classifiers": aug_classifiers}, "aug_pass")

# these cutoffs determine whether a transcript is pass/good/fail in addition to the classifiers
pass_cutoffs = {"coverage": 100.0, "percent_n": 1.0, "percent_coding_n": 0.2}
good_cutoffs = {"coverage": 95.0, "percent_n": 5.0, "percent_coding_n": 1.0}

tm_pass = dict_to_named_tuple(merge_dicts([{"classifiers": tm_coding_classifiers}, pass_cutoffs]), "tm_pass")
tm_good = dict_to_named_tuple(merge_dicts([{"classifiers": tm_coding_classifiers}, good_cutoffs]), "tm_good")

noncoding_pass = dict_to_named_tuple(merge_dicts([{"classifiers": noncoding_classifiers}, pass_cutoffs]),
                                     "noncoding_pass")
noncoding_good = dict_to_named_tuple(merge_dicts([{"classifiers": noncoding_classifiers}, good_cutoffs]),
                                     "noncoding_good")


# used for the plots
width = 9.0
height = 6.0
bar_width = 0.45
# paired_palette has two parallel color spectrums and black as the outgroup color
paired_palette = ["#df65b0", "#dd1c77", "#980043", "#a1dab4", "#41b6c4", "#2c7fb8", "#252525"]
# palette is the seaborn colorbind palette
palette = ["#0072b2", "#009e73", "#d55e00", "#cc79a7", "#f0e442", "#56b4e9"]


def refClassifiers(genome):
    base_query = "SELECT {} FROM details.'{}'"
    classifiers = classes_in_module(src.classifiers)
    classifiers = ",".join([x.__name__ for x in classifiers])
    query = base_query.format(classifiers, genome)
    return query


def allClassifiers(genome):
    base_query = "SELECT {} FROM details.'{}'"
    classifiers = classes_in_module(src.classifiers) + classes_in_module(src.alignment_classifiers)
    classifiers = ",".join([x.__name__ for x in classifiers])
    query = base_query.format(classifiers, genome)
    return query


def allAugustusClassifiers(genome):
    base_query = "SELECT {} FROM augustus_details.'{}'"
    classifiers = ",".join([x.__name__ for x in classes_in_module(src.augustus_classifiers)])
    query = base_query.format(classifiers, genome)
    return query


def potentiallyInterestingBiology(genome):
    base_query = ("SELECT details.'{0}'.InFrameStop,details.'{0}'.CodingMult3Insertions,"
                  "details.'{0}'.CodingMult3Deletions,details.'{0}'.Nonsynonymous,details.'{0}'.FrameShift FROM "
                  "details.'{0}' JOIN main.'{0}' USING ('AlignmentId') WHERE {1} ")
    equality = ["main.'{}'.{} = 0".format(genome, x) for x in tm_coding_classifiers]
    classifiers = " AND ".join(equality)
    query = base_query.format(genome, classifiers)
    return query


def assemblyErrors(genome, details=True):
    base_query = ("FROM details.'{0}' JOIN main.'{0}' USING "
                  "('AlignmentId') WHERE main.'{0}'.AlignmentPartialMap = 1 OR main.'{0}'.UnknownBases = 1 OR "
                  "main.'{0}'.UnknownGap = 1 OR main.'{0}'.ShortCds = 1 OR main.'{0}'.AlignmentAbutsUnknownBases = 1 "
                  "OR main.'{0}'.AlignmentAbutsRight = 1 OR main.'{0}'.AlignmentAbutsLeft = 1")
    details_selection = ("SELECT details.'{0}'.AlignmentPartialMap,details.'{0}'.UnknownBases,details.'{0}'.UnknownGap,"
                  "details.'{0}'.ShortCds,details.'{0}'.AlignmentAbutsUnknownBases,details.'{0}'.AlignmentAbutsRight,"
                  "details.'{0}'.AlignmentAbutsLeft ")
    classify_selection = ("SELECT main.'{0}'.AlignmentId ")
    added_query = details_selection + base_query if details else classify_selection + base_query
    query = added_query.format(genome)
    return query


def alignmentErrors(genome, details=True):
    base_query = ("FROM details.'{0}' JOIN main.'{0}' USING (AlignmentId) WHERE main.'{0}'.BadFrame = 1 OR "
                  "main.'{0}'.CdsGap = 1 OR main.'{0}'.CdsMult3Gap = 1 OR main.'{0}'.UtrGap = 1 OR "
                  "main.'{0}'.Paralogy = 1 OR main.'{0}'.HasOriginalIntrons = 1 OR main.'{0}'.StartOutOfFrame = 1")
    details_selection = ("SELECT details.'{0}'.BadFrame,details.'{0}'.CdsGap,details.'{0}'.CdsMult3Gap,"
                         "details.'{0}'.UtrGap,details.'{0}'.Paralogy,details.'{0}'.HasOriginalIntrons,"
                         "details.'{0}'.StartOutOfFrame ")
    classify_selection = "SELECT main.'{0}'.AlignmentId "
    added_query = details_selection + base_query if details else classify_selection + base_query
    query = added_query.format(genome)
    return query


def transMapEval(genome, coding=True, good=False):
    if coding is True:
        requirements = tm_pass if good is False else tm_good
    else:
        requirements = noncoding_pass if good is False else noncoding_good
    base_query = "SELECT main.'{0}'.AlignmentId FROM main.'{0}' JOIN attributes.'{0}' USING (AlignmentId) WHERE {1}"
    classifiers = ["main.'{}'.{} = 0".format(genome, x) for x in requirements.classifiers]
    classifiers = " AND ".join(classifiers)
    classifiers += " AND attributes.'{}'.AlignmentCoverage >= {}".format(genome, requirements.coverage)
    classifiers += " AND attributes.'{}'.PercentN < {}".format(genome, requirements.percent_n)
    if coding is True:
        classifiers += " AND attributes.'{}'.PercentCodingN < {}".format(genome, requirements.percent_coding_n)
    query = base_query.format(genome, classifiers)
    return query


def refEval(genome):
    """
    We only evaluate coding genes in the reference, otherwise we just don't have enough to say anything
    """
    base_query = "SELECT AlignmentId FROM main.'{0}' WHERE {1}"
    classifiers = ["main.'{}'.{} = 0".format(genome, x) for x in ref_pass.classifiers]
    classifiers = " AND ".join(classifiers)
    query = base_query.format(genome, classifiers)
    return query


def augustusEval(genome):
    base_query = "SELECT AlignmentId FROM augustus.'{0}' WHERE {1}"
    equality = ["{} = 0".format(x) for x in aug_pass.classifiers]
    classifiers = " AND ".join(equality)
    query = base_query.format(genome, classifiers)
    return query