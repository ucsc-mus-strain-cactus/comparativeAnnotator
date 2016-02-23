"""
This is the main driver script for comparativeAnnotator in transMap mode.
"""
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from pycbio.sys.dataOps import grouper
from pycbio.sys.introspection import classes_in_module
from pycbio.bio.transcripts import get_transcript_dict
from pycbio.bio.psl import get_alignment_dict
from pycbio.bio.bio import get_sequence_dict
from comparativeAnnotator.lib.name_conversions import remove_alignment_number
import comparativeAnnotator.reference_classifiers as ref_classifiers
from comparativeAnnotator.database_schema import construct_ref_tables

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
        ref_classifier_list = classes_in_module(ref_classifiers)
        # instantiate each classifier
        ref_classifier_fns = [x() for x in ref_classifier_list]
        ref_fasta = get_sequence_dict(self.args.ref_fasta)
        for aln_id, (a, t, aln) in self.chunk:
            r[aln_id] = {classify_fn.name: classify_fn(a, ref_fasta) for classify_fn in ref_classifier_fns}
            r[aln_id] = {classify_name: len(details) for classify_name, details in r[aln_id].iteritems()}
        df = pd.DataFrame.from_dict(r)
        df.to_sql(self.args.genome + '_Classify', self.engine, if_exists='append')
        df_details = pd.DataFrame.from_dict(r_details)
        df_details.to_sql(self.args.genome + '_Details', self.engine, if_exists='append')


def build_attributes_table(args, tx_dict, engine):
    """
    Produces the attributes table. Dumps the contents of the attributes.tsv file in addition to a few other metrics.
    """
    r = {}
    col_template = ['GeneId', 'GeneName', 'GeneType', 'TranscriptId', 'TranscriptType', 'NumberIntrons']
    with open(args.gencode_attributes) as f:
        for line in f:
            gene_id, gene_name, gene_type, transcript_id, transcript_type = line.split()
            num_introns = len(tx_dict[transcript_id].intron_intervals)
            d = dict(zip(*(col_template, (gene_id, gene_name, gene_type, transcript_id, transcript_type, num_introns))))
            r[transcript_id] = d
    df = pd.DataFrame.from_dict(r)
    df.to_sql(args.genome, engine, if_exists='append')


def chunk_ref_transcripts(target, args, engine, chunk_size=500):
    """
    Main loop for classification. Produces a classification job for chunk_size alignments.
    """
    def build_aln_dict(ref_dict, tx_dict, psl_dict):
        r = {}
        for aln_id, aln in psl_dict.iteritems():
            a = ref_dict[remove_alignment_number(aln_id)]
            t = tx_dict[aln_id]
            r[aln_id] = (a, t, aln)
        return r
    ref_dict = get_transcript_dict(args.annotation_gp)
    tx_dict = get_transcript_dict(args.target_gp)
    psl_dict = get_alignment_dict(args.psl)
    aln_dict = build_aln_dict(ref_dict, tx_dict, psl_dict)
    for chunk in grouper(aln_dict.iteritems(), chunk_size):
        target.addChildTarget(Classify(args, chunk, engine))
    target.addChildTargetFn(build_attributes_table, args=(args, psl_dict, tx_dict, engine))


def initialize_tables(target, args):
    ref_classify, ref_details, ref_attributes = construct_ref_tables(args.ref_genome)
    engine = create_engine('sqlite:///' + args.db)
    Base.metadata.create_all(engine)
    target.setFollowOnTargetFn(chunk_ref_transcripts, args=(args, engine))


def classify_startup(args):
    """
    Entry to start jobTree for classification.
    """
    target = Target.makeTargetFn(initialize_tables, args=(args,))
    failures = Stack(target).startJobTree(args)
    if failures != 0:
        raise Exception('Error: ' + str(failures) + ' jobs failed')
