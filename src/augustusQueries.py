import src.augustusClassifiers
from lib.general_lib import classesInModule

def augustusNotOk():
    """
    Defines augustus OK as not failing all Augustus classifiers.
    """
    classifyFields = detailsFields = [x.__name__ for x in classesInModule(src.augustusClassifiers)]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * (len(classifyFields))
    return detailsFields, classifyFields, classifyValues, classifyOperations