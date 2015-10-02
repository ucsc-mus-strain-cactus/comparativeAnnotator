import os
from itertools import izip_longest, product, izip
import sqlite3 as sql
from collections import defaultdict

import src.queries
import lib.sqlite_lib as sql_lib
import lib.psl_lib as psl_lib
import lib.sequence_lib as seq_lib
from lib.general_lib import functions_in_module
from src.abstractClassifier import AbstractClassifier
from jobTree.scriptTree.target import Target
from jobTree.src.bioio import logger, system

class BuildTracks(Target):