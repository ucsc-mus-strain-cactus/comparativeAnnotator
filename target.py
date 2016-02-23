"""
This is the main driver script for comparativeAnnotator in transMap mode.
"""
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from pycbio.sys.dataOps import grouper, merge_dicts
from pycbio.sys.introspection import classes_in_module
from pycbio.bio.transcripts import get_transcript_dict
from pycbio.bio.psl import get_alignment_dict
from pycbio.bio.bio import get_sequence_dict
from comparativeAnnotator.lib.name_conversions import remove_alignment_number
import comparativeAnnotator.reference_classifiers as ref_classifiers
import comparativeAnnotator.alignment_classifiers as aln_classifiers

__author__ = "Ian Fiddes"

Base = declarative_base()


class Classify(Target):
    """
    Main class for producing classifications.
    """
    def __init__(self, args, chunk, engine):
        Target.__init__(self)
        self.args = args
        self.chunk = chunk
        self.engine = engine

    def run(self):
        r = {}  # results dict contains integer classifications
        r_details = {}  # details dict contains BED entries
        # build list of classifier classes
        ref_classifier_list = classes_in_module(ref_classifiers)
        aln_classifier_list = classes_in_module(aln_classifiers)
        # instantiate, now they are callable
        aln_classifier_fns = [x() for x in ref_classifier_list]
        ref_classifier_fns = [x() for x in aln_classifier_list]
        # get pyfasta Fasta objects for each reference and target fastas
        ref_fasta = get_sequence_dict(self.args.ref_fasta)
        tgt_fasta = get_sequence_dict(self.args.fasta)
        for aln_id, (a, t, aln, ref_aln, paralogy_rec) in self.chunk:
            ref_results = {classify_fn.name: classify_fn(a, ref_fasta) for classify_fn in ref_classifier_fns}
            aln_results = {classify_fn.name: classify_fn(a, t, aln, ref_aln, ref_fasta, tgt_fasta) for classify_fn in
                           aln_classifier_fns}
            aln_results['Paralogy'] = paralogy_rec
            r_details[aln_id] = merge_dicts([ref_results, aln_results])
            r[aln_id] = {classify_name: len(details) for classify_name, details in r_details[aln_id].iteritems()}
        df = pd.DataFrame.from_dict(r)
        df.to_sql(self.args.genome, self.engine, if_exists='append')
        df_details = pd.DataFrame.from_dict(r_details)
        df_details.to_sql(self.args.genome + '_Details', self.engine, if_exists='append')


def chunk_tm_transcripts(target, args, engine, chunk_size=150):
    """
    Main loop for classification. Produces a classification job for chunk_size alignments.
    """
    def build_aln_dict(ref_dict, tx_dict, psl_dict, ref_psl_dict, paralogy_recs):
        r = {}
        for aln_id, aln in psl_dict.iteritems():
            a = ref_dict[remove_alignment_number(aln_id)]
            t = tx_dict[aln_id]
            ref_aln = ref_psl_dict[remove_alignment_number(aln_id)]
            c = paralogy_recs[aln_id]
            r[aln_id] = (a, t, aln, ref_aln, c)
        return r
    ref_dict = get_transcript_dict(args.annotation_gp)
    tx_dict = get_transcript_dict(args.target_gp)
    psl_dict = get_alignment_dict(args.psl)
    ref_psl_dict = get_alignment_dict(args.ref_psl)
    # we have to do paralogy separately, as it needs all the alignment names
    paralogy_recs = aln_classifiers.paralogy(tx_dict)
    aln_dict = build_aln_dict(ref_dict, tx_dict, psl_dict, ref_psl_dict, paralogy_recs)
    for chunk in grouper(aln_dict.iteritems(), chunk_size):
        target.addChildTarget(Classify(args, chunk, engine))


def classify_startup(args):
    """
    Entry to start jobTree for classification.
    """
    if args.mode == 'transMap':
        target = Target.makeTargetFn(construct_tm_tables, args=(args,))
    else:
        raise NotImplementedError
    failures = Stack(target).startJobTree(args)
    if failures != 0:
        raise Exception('Error: ' + str(failures) + ' jobs failed')
