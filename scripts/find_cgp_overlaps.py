"""
Find CGP overlaps, looking only at transcripts which are tagged as best
"""
import argparse
import sys
import os
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
import numpy as np
from collections import defaultdict, OrderedDict, Counter
from comparativeAnnotator.database_queries import get_row_dict, get_fail_pass_excel_ids, augustus_eval
from pycbio.sys.dataOps import merge_dicts
from pycbio.sys.mathOps import format_ratio
from pycbio.bio.transcripts import get_transcript_dict
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers, remove_augustus_alignment_number, \
    aln_id_is_augustus, aln_id_is_transmap
from comparativeAnnotator.database_queries import get_gene_transcript_map, get_transcript_gene_map, \
    get_transcript_biotype_map


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('cgp')
    parser.add_argument('tm')
