import os
import peewee
import pandas as pd
import cPickle as pickle
from collections import defaultdict

from jobTree.scriptTree.target import Target

from pycbio.sys.introspection import classes_in_module
from pycbio.bio.transcripts import get_transcript_dict
from pycbio.bio.psl import get_alignment_dict
from pycbio.bio.bio import get_sequence_dict
from pycbio.sys.dataOps import grouper
from pycbio.sys.fileOps import tmpFileGet, ensureDir
from pycbio.sys.sqliteOps import ExclusiveSqlConnection

from comparativeAnnotator.comp_lib.name_conversions import remove_alignment_number
from comparativeAnnotator.database_schema import tgt_tables, ref_tables, initialize_tables, fetch_database
import comparativeAnnotator.classifiers as classifiers
import comparativeAnnotator.alignment_classifiers as alignment_classifiers
#import comparativeAnnotator.augustus_classifiers as augustus_classifiers
import comparativeAnnotator.alignment_attributes as alignment_attributes

__author__ = "Ian Fiddes"


class AbstractClassify(Target):
    """
    Abstract class for producing classifications. Provides:
    1) A shared __init__
    2) A shared way to instantiate classifiers.
    3) A shared way to write to sql.
    """
    def __init__(self, args, chunk, tmp_dir):
        Target.__init__(self)
        self.args = args
        self.chunk = chunk
        self.tmp_dir = tmp_dir

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

    def write_classifications(self, r_details):
        """
        Write classifications to disk.
        """
        r_classify, r_details_string = self._munge_classify_dict(r_details)
        for d, table in [[r_classify, genome + '_Classify'], [r_details_string, genome + '_Details']]:
            self.write_to_disk(d, table)

    def write_attributes(self, attrs):
        df = self._munge_dict_to_dataframe(attrs)
        self.write_to_disk(df, tables.attrs)

    def write_to_disk(self, d, prefix=None):
        """
        Generic database writing function.
        """
        tmp_file = tmpFileGet(prefix=prefix, tmpDir=self.tmp_dir)
        with open(tmp_file, 'w') as outf:
            pickle.dump(d, outf)


class Classify(AbstractClassify):
    """
    Main class for single genome classifications.
    """
    def run(self):
        fasta = get_sequence_dict(self.args.ref_fasta)
        classifier_fns = self._instantiate_classifiers(classifiers)
        r_details = []
        for tx_id, a in self.chunk:
            rd = {'TranscriptId': tx_id}
            for classify_fn in classifier_fns:
                rd[classify_fn.name] = classify_fn(a, fasta)
            r_details.append(rd)
        r_classify, r_details_string = self._munge_classify_dict(r_details)
        for d, prefix in [[r_classify, 'classify'], [r_details_string, 'details']]:
            self.write_to_disk(d, prefix=prefix)


class AlignmentClassify(AbstractClassify):
    """
    Main class for alignment classifications.
    """
    def run(self):
        ref_fasta = get_sequence_dict(self.args.ref_fasta)
        tgt_fasta = get_sequence_dict(self.args.fasta)
        aln_classifier_fns = self._instantiate_classifiers(alignment_classifiers)
        classifier_fns = self._instantiate_classifiers(classifiers)
        r_details = []
        for aln_id, (a, t, aln, ref_aln, paralogy_rec) in self.chunk:
            rd = {'AlignmentId': aln_id, 'TranscriptId': a.name, 'Paralogy': paralogy_rec}
            for aln_classify_fn in aln_classifier_fns:
                rd[aln_classify_fn.name] = aln_classify_fn(a, t, aln, ref_aln, ref_fasta, tgt_fasta)
            for classify_fn in classifier_fns:
                rd[classify_fn.name] = classify_fn(t, tgt_fasta)
            r_details.append(rd)
        r_classify, r_details_string = self._munge_classify_dict(r_details)
        for d, prefix in [[r_classify, 'classify'], [r_details_string, 'details']]:
            self.write_to_disk(d, prefix=prefix)


class AlignmentAttributes(AbstractClassify):
    """
    Produces alignment attributes table.
    """
    def run(self):
        tables = tgt_tables(self.args.genome)
        ref_fasta = get_sequence_dict(self.args.ref_fasta)
        tgt_fasta = get_sequence_dict(self.args.fasta)
        attrs = self._instantiate_classifiers(alignment_attributes)
        r_attrs = {}
        for aln_id, (a, t, aln, ref_aln, paralogy_rec) in self.chunk:
            ra = {'TranscriptId': a.name}
            for attr_fn in attrs:
                ra[attr_fn.name] = attr_fn(a, t, aln, ref_aln, ref_fasta, tgt_fasta)
            r_attrs[aln_id] = ra
        self.write_to_disk(r_attrs)


def build_attributes_table(target, args, ref_dict, tmp_attrs):
    """
    Produces the reference attributes table.
    Dumps the contents of the attributes.tsv file in addition to reporting the number of introns and reference
    chromosome.
    """
    tables = ref_tables(args.ref_genome)
    attr_table_name = tables.attrs.__name__
    df = pd.read_table(args.gencode_attributes, sep='\t', index_col=3, header=0)
    d = []
    for ens_id, a in ref_dict.iteritems():
        row = {'TranscriptId': ens_id}
        row['RefChrom'] = a.chromosome
        row['NumberIntrons'] = len(a.intron_intervals)
        row.update(df.loc[ens_id].to_dict())
        d.append(row)
    tmp_file = tmpFileGet(tmpDir=tmp_attrs)
    with open(tmp_file, 'w') as outf:
        pickle.dump(d, outf)


def write_to_db(target, args, tmp_classify, tmp_attrs, primary_key):
    """
    Wrapper to find pickled objects and write them to the database.
    """
    def peewee_write(table, db, files):
        """
        We have to write in chunks to avoid OperationalError: too many SQL variables
        the limit is 999, and it is really annoying to change.
        TODO: you could consider producing many small database files in parallel and then merging using CLI
        """
        for f in files:
            data = pickle.load(open(f))
            with db.atomic():
                for idx in range(0, len(data), 50):
                    table.insert_many(data[idx:idx + 50]).execute()
    if args.mode == 'reference':
        tables = ref_tables(args.ref_genome)
    elif args.mode == 'transMap':
        tables = tgt_tables(args.genome)
    initialize_tables(tables, args.db)
    db = fetch_database()
    attr_pickle_files = [os.path.join(tmp_attrs, x) for x in os.listdir(tmp_attrs)]
    peewee_write(tables.attrs, db, attr_pickle_files)
    details_pickle_files = []
    classify_pickle_files = []
    for f in os.listdir(tmp_classify):
        if 'details' in f:
            details_pickle_files.append(os.path.join(tmp_classify, f))
        elif 'classify' in f:
            classify_pickle_files.append(os.path.join(tmp_classify, f))
        else:
            raise Exception("Ian is a bad programmer")
    assert len(details_pickle_files) == len(classify_pickle_files)
    for table, files in [[tables.classify, classify_pickle_files], [tables.details, details_pickle_files]]:
        peewee_write(table, db, files)


def construct_tmp_dirs(target):
    """
    Shared tmp_dir construction for reference or tm
    """
    tmp_dir = target.getGlobalTempDir()
    tmp_classify = os.path.join(tmp_dir, 'classify')
    tmp_attrs = os.path.join(tmp_dir, 'attrs')
    ensureDir(tmp_classify)
    ensureDir(tmp_attrs)
    return tmp_classify, tmp_attrs


def run_ref_classifiers(target, args, primary_key='TranscriptId', chunk_size=2000):
    """
    Main loop for classification. Produces a classification job for chunk_size alignments.
    """
    tmp_classify, tmp_attrs = construct_tmp_dirs(target)
    ref_dict = get_transcript_dict(args.annotation_gp)
    for chunk in grouper(ref_dict.iteritems(), chunk_size):
        target.addChildTarget(Classify(args, chunk, tmp_classify))
    target.addChildTargetFn(build_attributes_table, args=(args, ref_dict, tmp_attrs))
    target.setFollowOnTargetFn(write_to_db, args=(args, tmp_classify, tmp_attrs, primary_key))


def run_tm_classifiers(target, args, primary_key='AlignmentId', chunk_size=150):
    """
    Main loop for classification. Produces a classification job for chunk_size alignments.
    """
    tmp_classify, tmp_attrs = construct_tmp_dirs(target)
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
        target.addChildTarget(AlignmentClassify(args, chunk, tmp_classify))
        target.addChildTarget(AlignmentAttributes(args, chunk, tmp_attrs))
    target.setFollowOnTargetFn(write_to_db, args=(args, tmp_classify, tmp_attrs, primary_key))


def run_aug_classifiers(target, args, primary_key='AlignmentId', chunk_size=150):
    """
    Main loop for augustus classification. Produces a classification job for chunk_size alignments.
    """
    tmp_classify, tmp_attrs = construct_tmp_dirs(target)
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
        target.addChildTarget(AlignmentClassify(args, chunk, tmp_classify))
        target.addChildTarget(AlignmentAttributes(args, chunk, tmp_attrs))
    target.setFollowOnTargetFn(write_to_db, args=(args, tmp_classify, tmp_attrs, primary_key))