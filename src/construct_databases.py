"""
Merges all of the pickled dicts into a sqlite table
"""
import os
import cPickle as pickle
import lib.sql_lib as sql_Lib
from jobTree.scriptTree.target import Target


class ConstructDatabases(Target):
    def __init__(self, genome, out_classify, out_details, out_attributes):
        Target.__init__(self)
        self.genome = genome
        self.out_classify = out_classify
        self.out_details = out_details
        self.out_attributes = out_attributes

    def load_data(self):
        for t in ["Classify", "Details", "Attribute"]:
            data_dict = {}
            for f in os.listdir(self.getGlobalTempDir()):
                if t in f:
                    col = f.split(t)[1]
                    assert col not in data_dict
                    data_dict[col] = pickle.load(f)
            yield data_dict

    def run(self):
        for d in self.load_data():
            sql_lib.write_dict(d)