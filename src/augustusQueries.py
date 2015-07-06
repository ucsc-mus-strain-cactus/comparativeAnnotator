import src.classifiers
from lib.general_lib import classesInModule

def augustusOk():
    """
    Defines augustus OK as not failing all Augustus classifiers.
    """
    classifyFields = detailsFields = [x.__name__ for x in classesInModule(src.augustusClassifiers)]
    classifyOperations = ["AND"] * (len(classifyFields) - 1)
    classifyValues = [1] * (len(classifyFields) - 1)
    return detailsFields, classifyFields, classifyValues, classifyOperations