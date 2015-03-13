
def mutations():
    """
    Hunts for likely real mutations by excluding assembly and alignment errors.
    Any transcript which has errors in the classify fields specified will not have the details shown.
    """
    detailsFields = ["CodingInsertions", "CodingDeletions", "CodingMult3Insertions", "CodingMult3Deletions",
                     "CdsNonCanonSplice", "UtrNonCanonSplice", "CdsUnknownSplice", "UtrUnknownSplice", "CdsMult3Gap", "InFrameStop",
                     "Synonymous", "Nonsynonymous", "Frameshift"]
    classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases", "ScaffoldGap", "NoCds", 
                     "AlignmentPartialMap", "BadFrame", "CdsGap", "UtrGap"]
    classifyOperations = ["AND"] * (len(classifyFields) - 1)
    classifyValues = [0] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations

def inFrameStop():
    detailsFields = ["InFrameStop"]
    classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases", "ScaffoldGap", "NoCds", 
                     "AlignmentPartialMap", "BadFrame", "CdsGap", "UtrGap"]
    classifyOperations = ["AND"] * (len(classifyFields) - 1)
    classifyValues = [0] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations

def interestingBiology():
    detailsFields = ["InFrameStop", "CodingMult3Insertions", "CodingMult3Deletions", "CdsMult3Gap", "Frameshift"]
    classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases", "ScaffoldGap", "NoCds", 
                     "AlignmentPartialMap", "BadFrame", "CdsGap", "UtrGap"]
    classifyOperations = ["AND"] * (len(classifyFields) - 1)
    classifyValues = [0] * len(classifyFields)  
    return detailsFields, classifyFields, classifyValues, classifyOperations

def assemblyErrors():
    """
    Looks for assembly errors. Reports transcripts with assembly errors.
    """
    detailsFields = classifyFields = ["AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "UnknownBases", "ScaffoldGap", "NoCds"]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations

def alignmentErrors():
    """
    Looks for alignment errors. Reports details for all fields that are likely alignment errors.
    """
    classifyFields = detailsFields = ["BadFrame", "CdsGap", "UtrGap"]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * len(classifyFields)
    return detailsFields, classifyFields, classifyValues, classifyOperations
