"""
This file declares configuration variables used for plotting.
"""

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
aug_ok_fields = ['AugustusParalogy', 'AugustusExonGain', 'AugustusExonLoss', 'AugustusNotSameStrand',
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
