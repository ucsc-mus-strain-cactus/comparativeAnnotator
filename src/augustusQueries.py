import src.augustusClassifiers
from lib.general_lib import classesInModule

def augustusOk():
    """
    Defines augustus OK as not failing all Augustus classifiers.
    """
    classifyFields = detailsFields = [x.__name__ for x in classesInModule(src.augustusClassifiers)]
    classifyOperations = ["AND"] * (len(classifyFields) - 1)
    classifyValues = [0] * (len(classifyFields))
    return detailsFields, classifyFields, classifyValues, classifyOperations