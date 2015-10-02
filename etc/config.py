"""
This file declares configuration variables used for plotting.
"""

import src.classifiers
import src.augustus_classifiers
from lib.general_lib import classes_in_module

__author__ = "Ian Fiddes"


# Below are global variables shared by plotting scripts

# this genetic distance is from Joel. There are functions to derive it, but this makes it consistent across plots.
hard_coded_genome_order = ['C57B6NJ', 'NZOHlLtJ', '129S1', 'FVBNJ', 'NODShiLtJ', 'LPJ', 'AJ', 'AKRJ', 'BALBcJ', 'DBA2J',
                            'C3HHeJ', 'CBAJ', 'WSBEiJ', 'CASTEiJ', 'PWKPhJ', 'SPRETEiJ', 'CAROLIEiJ', 'PAHARIEiJ']

# these classifiers define OK for coding transcripts
tm_coding_classifiers = ["CodingInsertions", "CodingDeletions",
                         "AlignmentPartialMap", "BadFrame", "BeginStart", "UnknownBases", "AlignmentAbutsUnknownBases",
                         "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice",
                         "EndStop", "InFrameStop", "ShortCds", "StartOutOfFrame", "FrameShift",
                         "AlignmentAbutsRight", "AlignmentAbutsLeft"]

# these classifiers define OK for non-coding transcripts
tm_noncoding_classifiers = ["AlignmentPartialMap", "UtrUnknownSplice", "UtrGap", "UnknownGap", "UnknownBases",
                            "AlignmentAbutsUnknownBases"]

# these classifiers define OK for Augustus transcripts
aug_ok_classifiers = ['AugustusParalogy', 'AugustusExonGain', 'AugustusExonLoss', 'AugustusNotSameStrand',
                      'AugustusNotSameStartStop', 'AugustusNotSimilarTerminalExonBoundaries',
                      'AugustusNotSimilarInternalExonBoundaries']

# used for the plots
width = 9.0
height = 6.0
bar_width = 0.45
# paired_palette has two parallel color spectrums and black as the outgroup color
paired_palette = ["#df65b0", "#dd1c77", "#980043", "#a1dab4", "#41b6c4", "#2c7fb8", "#252525"]
# palette is the seaborn colorbind palette
palette = ["#0072b2", "#009e73", "#d55e00", "#cc79a7", "#f0e442", "#56b4e9"]


def allClassifiers(genome):
    base_query = "SELECT {} FROM details.'{}'"
    classifiers = ",".join([x.__name__ for x in classes_in_module(src.classifiers)])
    query = base_query.format(classifiers, genome)
    return query


def allAugustusClassifiers(genome):
    base_query = "SELECT {} FROM details.'{}'"
    classifiers = ",".join([x.__name__ for x in classes_in_module(src.augustus_classifiers)])
    query = base_query.format(classifiers, genome)
    return query


def potentiallyInterestingBiology(genome):
    base_query = ("SELECT details.InFrameStop,details.CodingMult3Insertions,details.CodingMult3Deletions,"
                  "details.Nonsynonymous,details.FrameShift FROM details.'{0}' JOIN main.'{0}' USING ('AlignmentId') "
                  "WHERE {1} ")
    equality = ["main.{} = 0".format(x) for x in tm_coding_classifiers]
    classifiers = " AND ".join(equality)
    query = base_query.format(genome, classifiers)
    return query


def assemblyErrors(genome):
    base_query = ("SELECT details.AlignmentPartialMap,details.UnknownBases,details.UnknownGap,details.ShortCds,"
                  "details.AlignmentAbutsUnknownBases,details.AlignmentAbutsRight,details.AlignmentAbutsLeft "
                  "FROM details.'{0}' JOIN main.'{0}' USING ('AlignmentId') WHERE main.AlignmentPartialMap = 1 "
                  "OR main.UnknownBases = 1 OR main.UnknownGap = 1 OR main.ShortCds = 1 OR "
                  "main.AlignmentAbutsUnknownBases = 1 OR main.AlignmentAbutsRight = 1 OR main.AlignmentAbutsLeft = 1")
    query = base_query.format(genome)
    return query


def alignmentErrors(genome):
    base_query = ("SELECT details.BadFrame,details.CdsGap,details.CdsMult3Gap,details.UtrGap,details.Paralogy,"
                  "details.HasOriginalIntrons,details.StartOutOfFrame FROM details.'{0}' JOIN main.'{0}' USING "
                  "(AlignmentId) WHERE main.BadFrame = 1 OR main.CdsGap = 1 OR main.CdsMult3Gap = 1 OR "
                  "main.UtrGap = 1 OR main.Paralogy = 1 OR main.HasOriginalIntrons = 1")
    query = base_query.format(genome)
    return query


def transMapOk(genome):
    base_query = "SELECT AlignmentId FROM main.'{}' WHERE {}"
    equality = ["{} = 0".format(x) for x in tm_coding_classifiers]
    classifiers = " AND ".join(equality)
    query = base_query.format(genome, classifiers)
    return query


def augustusOk(genome):
    base_query = "SELECT AlignmentId FROM main.'{}' WHERE {}"
    equality = ["{} = 0".format(x) for x in aug_ok_classifiers]
    classifiers = " AND ".join(equality)
    query = base_query.format(genome, classifiers)
    return query