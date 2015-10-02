import src.augustus_classifiers
from lib.general_lib import classes_in_module


def augustus_not_ok():
    """
    Defines augustus OK as not failing all Augustus classifiers.
    """
    classify_fields = details_fields = [x.__name__ for x in classes_in_module(src.augustus_classifiers)]
    classify_operations = ["OR"] * (len(classify_fields) - 1)
    classify_values = [1] * (len(classify_fields))
    return details_fields, classify_fields, classify_values, classify_operations