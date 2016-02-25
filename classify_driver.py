import os
import peewee
import pandas as pd
from collections import defaultdict

from jobTree.scriptTree.target import Target

from pycbio.sys.introspection import classes_in_module
from pycbio.bio.transcripts import get_transcript_dict
from pycbio.bio.psl import get_alignment_dict
from pycbio.bio.bio import get_sequence_dict
from pycbio.sys.dataOps import grouper

from comparativeAnnotator.comp_lib.name_conversions import remove_alignment_number
from comparativeAnnotator.database_schema import tgt_tables, ref_tables, initialize_tables
import comparativeAnnotator.classifiers as classifiers
import comparativeAnnotator.alignment_classifiers as alignment_classifiers
#import comparativeAnnotator.augustus_classifiers as augustus_classifiers
import comparativeAnnotator.alignment_attributes as alignment_attributes

__author__ = "Ian Fiddes"

recs_per_upsert = 25
# TODO: why is this value so low? Anything bigger leads to errors.


class AbstractClassify(Target):
    """
    Abstract class for producing classifications. Provides:
    1) A shared __init__
    2) A shared way to instantiate classifiers.
    3) A shared way to write to sql.
    """
    def __init__(self, args, chunk):
        Target.__init__(self)
        self.args = args
        self.chunk = chunk

    def _munge_classify_dict(self, r_details):
        """
        Creates a count from a details dict as well as converting lists to strings.
        """
        r_classify = []
        r_details_string = []
        for rd in r_details:
            rcs = {}
            rds = {}
            for classify_name, details in rd.iteritems():
                if isinstance(details, list):
                    rcs[classify_name] = len(details)
                    rds[classify_name] = '\n'.join(['\t'.join(map(str, bed)) for bed in details if len(bed) > 0])
                else:
                    assert isinstance(details, str)
                    rds[classify_name] = rcs[classify_name] = details
            r_details_string.append(rds)
            r_classify.append(rcs)
        return r_classify, r_details_string

    def _instantiate_classifiers(self, classifiers):
        """
        Instantiate each classifier class, turning them into function-like objects
        """
        classifier_list = classes_in_module(classifiers)
        return [c() for c in classifier_list]

    def write_classifications(self, r_details, genome, tables, primary_key):
        """
        Write classifications to the database.
        """
        r_classify, r_details_string = self._munge_classify_dict(r_details)
        for d, table in [[r_classify, tables.classify], [r_details_string, tables.details]]:
            self.write_to_db(d, table, primary_key)

    def write_attributes(self, attrs, genome, tables, primary_key):
        df = self._munge_dict_to_dataframe(attrs, primary_key)
        self.write_to_db(df, tables.attrs, primary_key)

    def write_to_db(self, d, table, primary_key):
        """
        Generic database writing function.
        """
        table._meta.database.init(self.args.db)
        with table._meta.database.atomic() as txn:
            for chunk in grouper(d, recs_per_upsert):
                _ = table.insert_many(chunk).upsert().execute()


class Classify(AbstractClassify):
    """
    Main class for single genome classifications.
    """
    def run(self, primary_key='TranscriptId'):
        tables = ref_tables(self.args.ref_genome)
        fasta = get_sequence_dict(self.args.ref_fasta)
        classifier_fns = self._instantiate_classifiers(classifiers)
        r_details = []
        for tx_id, a in self.chunk:
            rd = {primary_key: tx_id}
            for classify_fn in classifier_fns:
                rd[classify_fn.name] = classify_fn(a, fasta)
            r_details.append(rd)
        self.write_classifications(r_details, self.args.ref_genome, tables, primary_key)


class AlignmentClassify(AbstractClassify):
    """
    Main class for alignment classifications.
    """
    def run(self, primary_key='AlignmentId'):
        tables = tgt_tables(self.args.genome)
        ref_fasta = get_sequence_dict(self.args.ref_fasta)
        tgt_fasta = get_sequence_dict(self.args.fasta)
        aln_classifier_fns = self._instantiate_classifiers(alignment_classifiers)
        classifier_fns = self._instantiate_classifiers(classifiers)
        r_details = []
        for aln_id, (a, t, aln, ref_aln, paralogy_rec) in self.chunk:
            rd = {primary_key: aln_id, 'TranscriptId': a.name, 'Paralogy': paralogy_rec}
            for aln_classify_fn in aln_classifier_fns:
                rd[aln_classify_fn.name] = aln_classify_fn(a, t, aln, ref_aln, ref_fasta, tgt_fasta)
            for classify_fn in classifier_fns:
                rd[classify_fn.name] = classify_fn(t, tgt_fasta)
            r_details.append(rd)
        self.write_classifications(r_details, self.args.genome, tables, primary_key)


class AlignmentAttributes(AbstractClassify):
    """
    Produces alignment attributes table.
    """
    def run(self, primary_key='AlignmentId'):
        tables = tgt_tables(self.args.genome)
        ref_fasta = get_sequence_dict(self.args.ref_fasta)
        tgt_fasta = get_sequence_dict(self.args.fasta)
        attrs = self._instantiate_classifiers(alignment_attributes)
        r_attrs = {}
        for aln_id, (a, t, aln, ref_aln, paralogy_rec) in self.chunk:
            ra = {primary_key: aln_id, 'TranscriptId': a.name}
            for attr_fn in attrs:
                ra[attr_fn.name] = attr_fn(a, t, aln, ref_aln, ref_fasta, tgt_fasta)
            r_attrs[aln_id] = ra
        self.write_attributes(r_attrs, self.args.genome, tables, primary_key)


def build_attributes_table(target, args, ref_dict):
    """
    Produces the reference attributes table.
    Dumps the contents of the attributes.tsv file in addition to reporting the number of introns and reference
    chromosome.
    """
    tables = ref_tables(args.ref_genome)
    df = pd.read_table(args.gencode_attributes, sep='\t', index_col=3, header=0)
    assert set(df.columns) == set(['GeneId', 'GeneName', 'GeneType', 'TranscriptId', 'TranscriptType'])
    more_cols = defaultdict(dict)
    for ens_id, a in ref_dict.iteritems():
        more_cols['RefChrom'][ens_id] = a.chromosome
        more_cols['NumberIntrons'][ens_id] = len(a.intron_intervals)
    more_cols_df = pd.DataFrame.from_dict(more_cols)
    df_combined = pd.merge(df, more_cols_df, left_index=True, right_index=True)
    # need a list of dicts for peewee
    df_back_to_dict = [x[1].to_dict() for x in df.iterrows()]
    for chunk in grouper(df_back_to_dict, recs_per_upsert):
        _ = tables.attrs.insert_many(chunk).upsert().execute()


def run_ref_classifiers(target, args, chunk_size=2000):
    """
    Main loop for classification. Produces a classification job for chunk_size alignments.
    """
    tables = ref_tables(args.ref_genome)
    initialize_tables(tables, args.db)
    ref_dict = get_transcript_dict(args.annotation_gp)
    for chunk in grouper(ref_dict.iteritems(), chunk_size):
        target.addChildTarget(Classify(args, chunk))
    target.addChildTargetFn(build_attributes_table, args=(args, ref_dict))


def run_tm_classifiers(target, args, chunk_size=150):
    """
    Main loop for classification. Produces a classification job for chunk_size alignments.
    """
    tables = tgt_tables(args.genome)
    initialize_tables(tables, args.db)
    def build_aln_dict(ref_dict, tx_dict, psl_dict, ref_psl_dict, paralogy_recs):
        """merge different data dicts"""
        r = {}
        for aln_id, aln in psl_dict.iteritems():
            if aln_id not in tx_dict:
                # not all alignments have transcripts
                continue
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
    paralogy_recs = alignment_classifiers.paralogy(tx_dict)
    aln_dict = build_aln_dict(ref_dict, tx_dict, psl_dict, ref_psl_dict, paralogy_recs)
    for chunk in grouper(aln_dict.iteritems(), chunk_size):
        target.addChildTarget(AlignmentClassify(args, chunk))
        target.addChildTarget(AlignmentAttributes(args, chunk))
