import src.augustusClassifiers
from lib.general_lib import classes_in_module

def augustusNotOk():
    """
    Defines augustus OK as not failing all Augustus classifiers.
    """
    classifyFields = detailsFields = [x.__name__ for x in classes_in_module(src.augustusClassifiers)]
    classifyOperations = ["OR"] * (len(classifyFields) - 1)
    classifyValues = [1] * (len(classifyFields))
    return detailsFields, classifyFields, classifyValues, classifyOperations