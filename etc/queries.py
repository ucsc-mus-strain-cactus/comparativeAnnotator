import src.classifiers
from lib.general_lib import classes_in_module


def inFrameStop():
    detailsFields = ["InFrameStop"]
    classifyFields = ["AlignmentPartialMap", "UnknownBases", "ShortCds", "BadFrame", "CdsGap", "UtrGap"]
    classifyOperations = ["AND"] * (len(classifyFields) - 1)
    classifyValues = [0] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations


def allProblems():
    detailsFields = classifyFields = [x.__name__ for x in classes_in_module(src.classifiers)]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations


def interestingBiology():
    detailsFields = ["InFrameStop", "CodingMult3Insertions", "CodingMult3Deletions", "StartOutOfFrame", "Nonsynonymous", "FrameShift"]
    excludedClassifyFields = ["AlignmentPartialMap", "UnknownBases", "ShortCds", "AlignmentPartialMap", 
                              "BadFrame", "CdsGap", "UtrGap", "UnknownGap", "HasOriginalIntrons"]
    excludedClassifyValues = [0] * len(excludedClassifyFields)
    excludedClassifyOperations = ["AND"] * len(excludedClassifyFields)
    includedClassifyFields = detailsFields
    includedClassifyValues = [1] * len(includedClassifyFields)
    includedClassifyOperations = ["OR"] * (len(includedClassifyFields) - 1)
    classifyFields = includedClassifyFields + excludedClassifyFields
    classifyValues = includedClassifyValues + excludedClassifyValues
    classifyOperations = includedClassifyOperations + excludedClassifyOperations
    return detailsFields, classifyFields, classifyValues, classifyOperations


def assemblyErrors():
    """
    Looks for assembly errors. Reports transcripts with assembly errors.
    """
    detailsFields = classifyFields = ["AlignmentPartialMap", "UnknownBases", "UnknownGap", "ShortCds",
                                      "AlignmentAbutsUnknownBases", "AlignmentAbutsRight", "AlignmentAbutsLeft"]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations


def alignmentErrors():
    """
    Looks for alignment errors. Reports details for all fields that are likely alignment errors.
    """
    classifyFields = detailsFields = ["BadFrame", "CdsGap", "CdsMult3Gap", "UtrGap", "Paralogy", "HasOriginalIntrons"]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations
