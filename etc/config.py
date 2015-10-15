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
ref_coding_classifiers = ["BadFrame", "BeginStart", "EndStop", "CdsGap", "CdsUnknownSplice", "UtrUnknownSplice",
                          "StartOutOfFrame", "SpliceContainsUnknownBases", "InFrameStop", "ShortCds"]

# these classifiers define Pass for coding transcripts
tm_pass_classifiers = ["BadFrame", "BeginStart", "EndStop", "CdsGap", "CdsUnknownSplice", "UtrUnknownSplice",
                       "StartOutOfFrame", "InFrameStop", "ShortCds", "CodingInsertions", "CodingDeletions",
                       "FrameShift", "HasOriginalStart", "HasOriginalStop", "HasOriginalIntrons"]

# these classifiers define Good for coding transcripts
# the difference: Can have an incomplete CDS, but that incomplete CDS should remain in frame. UtrUnknownSplice is also
# allowed.
tm_good_classifiers = ["CdsUnknownSplice", "FrameShift", "CodingInsertions", "CodingDeletions", "HasOriginalIntrons"]

# these classifiers define Pass/Good for non-coding transcripts
noncoding_pass_classifiers = ref_noncoding_classifiers = ["UtrUnknownSplice"]
noncoding_good_classifiers = ["UtrUnknownSplice", "UtrGap"]

# these classifiers define Pass for Augustus transcripts
aug_classifiers = ['AugustusParalogy', 'AugustusExonGain', 'AugustusExonLoss', 'AugustusNotSameStrand',
                   'AugustusNotSameStartStop', 'AugustusNotSimilarTerminalExonBoundaries',
                   'AugustusNotSimilarInternalExonBoundaries']

ref_pass = dict_to_named_tuple({"classifiers": ref_coding_classifiers}, "ref_pass")
aug_pass = dict_to_named_tuple({"classifiers": aug_classifiers}, "aug_pass")

# these cutoffs determine whether a transcript is pass/good/fail in addition to the classifiers
pass_cutoffs = {"coverage": 100.0, "percent_n": 1.0, "percent_coding_n": 0.2}
good_cutoffs = {"coverage": 95.0, "percent_n": 5.0, "percent_coding_n": 1.0}

# we use the ref classifiers to not unfairly impinge problematic source transcripts
ref_coding_classifiers = {"ref_classifiers": ref_coding_classifiers}
ref_noncoding_classifiers = {"ref_classifiers": ref_noncoding_classifiers}

tm_pass = dict_to_named_tuple(merge_dicts([{"classifiers": tm_pass_classifiers}, pass_cutoffs, ref_coding_classifiers]),
                              "tm_pass")
tm_good = dict_to_named_tuple(merge_dicts([{"classifiers": tm_good_classifiers}, good_cutoffs, ref_coding_classifiers]),
                              "tm_good")

noncoding_pass = dict_to_named_tuple(merge_dicts([{"classifiers": noncoding_pass_classifiers},
                                                  pass_cutoffs, ref_noncoding_classifiers]), "noncoding_pass")
noncoding_good = dict_to_named_tuple(merge_dicts([{"classifiers": noncoding_good_classifiers},
                                                  good_cutoffs, ref_noncoding_classifiers]), "noncoding_good")


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
    equality = ["main.'{}'.{} = 0".format(genome, x) for x in tm_pass_classifiers]
    classifiers = " AND ".join(equality)
    query = base_query.format(genome, classifiers)
    return query


def assemblyErrors(genome, biotype, details=True):
    base_query = ("FROM attributes.'{0}' JOIN details.'{0}' USING ('AlignmentId') JOIN main.'{0}' USING "
                  "('AlignmentId') WHERE main.'{0}'.AlignmentPartialMap = 1 OR main.'{0}'.UnknownBases = 1 OR "
                  "main.'{0}'.UnknownGap = 1 OR main.'{0}'.ShortCds = 1 OR main.'{0}'.AlnAbutsUnknownBases = 1 "
                  "OR main.'{0}'.AlnExtendsOffContig = 1")
    details_selection = ("SELECT details.'{0}'.AlignmentPartialMap,details.'{0}'.UnknownBases,details.'{0}'.UnknownGap,"
                         "details.'{0}'.ShortCds,details.'{0}'.AlnAbutsUnknownBases,details.'{0}'.AlnExtendsOffContig")
    classify_selection = "SELECT main.'{0}'.AlignmentId "
    added_query = details_selection + base_query if details else classify_selection + base_query
    query = added_query.format(genome)
    query += " AND attributes.'{}'.TranscriptType = '{}'".format(genome, biotype)
    return query


def alignmentErrors(genome, biotype, details=True):
    base_query = ("FROM attributes.'{0}' JOIN details.'{0}' USING ('AlignmentId') JOIN main.'{0}' USING "
                  "('AlignmentId') WHERE main.'{0}'.BadFrame = 1 OR "
                  "main.'{0}'.CdsGap = 1 OR main.'{0}'.CdsMult3Gap = 1 OR main.'{0}'.UtrGap = 1 OR "
                  "main.'{0}'.Paralogy = 1 OR main.'{0}'.HasOriginalIntrons = 1 OR main.'{0}'.StartOutOfFrame = 1")
    details_selection = ("SELECT details.'{0}'.BadFrame,details.'{0}'.CdsGap,details.'{0}'.CdsMult3Gap,"
                         "details.'{0}'.UtrGap,details.'{0}'.Paralogy,details.'{0}'.HasOriginalIntrons,"
                         "details.'{0}'.StartOutOfFrame ")
    classify_selection = "SELECT main.'{0}'.AlignmentId "
    added_query = details_selection + base_query if details else classify_selection + base_query
    query = added_query.format(genome)
    query += " AND attributes.'{}'.TranscriptType = '{}'".format(genome, biotype)
    return query


def transMapEval(ref_genome, genome, biotype, good=False):
    if biotype == "protein_coding":
        requirements = tm_pass if good is False else tm_good
    else:
        requirements = noncoding_pass if good is False else noncoding_good
    base_query = ("SELECT AlignmentId FROM attributes.'{ref_genome}' JOIN main.'{ref_genome}' USING (TranscriptId) "
                  "JOIN attributes.'{genome}' USING (TranscriptId) JOIN main.'{genome}' USING "
                  "(AlignmentId) WHERE {classifiers}")
    classifiers = []
    for c in requirements.classifiers:
        if c in requirements.ref_classifiers:
            l = "NOT (main.'{ref_genome}'.{classifier} = 0 AND main.'{genome}'.{classifier} = 1)"
            l = l.format(ref_genome=ref_genome, classifier=c, genome=genome)
        else:
            l = "main.'{}'.{} = 0".format(genome, c)
        classifiers.append(l)
    classifiers = " AND ".join(classifiers)
    classifiers += " AND attributes.'{}'.AlignmentCoverage >= {}".format(genome, requirements.coverage)
    classifiers += " AND attributes.'{}'.PercentUnknownBases < {}".format(genome, requirements.percent_n)
    if biotype == "protein_coding":
        classifiers += " AND attributes.'{}'.PercentUnknownCodingBases < {}".format(genome,
                                                                                    requirements.percent_coding_n)
    classifiers += " AND attributes.'{}'.TranscriptType = '{}'".format(genome, biotype)
    query = base_query.format(genome=genome, ref_genome=ref_genome, classifiers=classifiers)
    return query


def refEval(genome):
    base_query = "SELECT TranscriptId FROM main.'{0}' WHERE {1}"
    classifiers = ["{} = 0".format(x) for x in ref_pass.classifiers]
    classifiers = " AND ".join(classifiers)
    query = base_query.format(genome, classifiers)
    return query


def augustusEval(genome):
    base_query = "SELECT AlignmentId FROM augustus.'{0}' WHERE {1}"
    equality = ["{} = 0".format(x) for x in aug_pass.classifiers]
    classifiers = " AND ".join(equality)
    query = base_query.format(genome, classifiers)
    return query