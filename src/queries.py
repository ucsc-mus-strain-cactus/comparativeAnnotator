import src.classifiers
from lib.general_lib import classesInModule

def inFrameStop():
    detailsFields = ["InFrameStop"]
    classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases", "ScaffoldGap", "ShortCds", 
                     "AlignmentPartialMap", "BadFrame", "CdsGap", "UtrGap"]
    classifyOperations = ["AND"] * (len(classifyFields) - 1)
    classifyValues = [0] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations

def everything():
    detailsFields = classifyFields = [x.__name__ for x in classesInModule(src.classifiers)]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations

def interestingBiology():
    detailsFields = ["InFrameStop", "CodingMult3Insertions", "CodingMult3Deletions", "StartOutOfFrame", "Nonsynonymous", "FrameShift"]
    excludedClassifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases", "ScaffoldGap", "ShortCds", 
                              "AlignmentPartialMap", "BadFrame", "CdsGap", "UtrGap", "UnknownGap"]
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
    detailsFields = classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases", "ScaffoldGap", "ShortCds", "UnknownGap"]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations

def alignmentErrors():
    """
    Looks for alignment errors. Reports details for all fields that are likely alignment errors.
    """
    classifyFields = detailsFields = ["BadFrame", "CdsGap", "CdsMult3Gap", "UtrGap", "Paralogy", "StartOutOfFrame"]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations
