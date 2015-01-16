import os

from jobTree.scriptTree.target import Target
from lib.general_lib import classesInModule
import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib
import lib.sqlite_lib as sql_lib
import src.classifiers


class ClassifierMaster(Target):
    def __init__(self, genome, alnPsl, seqTwoBit, refSeqTwoBit, annotationBed,   
                gencodeAttributeMap, geneCheckBed, outDir, refGenome, primaryKey):
        #initialize the Target
        Target.__init__(self)

        #store the output database location
        self.db = os.path.join(outDir, genome + ".db")

        #primary key this will be keyed on (alignmentID usually)
        self.primaryKey = primaryKey
        self.genome = genome
        self.refGenome = refGenome
        self.alnPsl = alnPsl
        self.seqTwoBit = seqTwoBit
        self.refSeqTwoBit = refSeqTwoBit
        self.gencodeAttributeMap = gencodeAttributeMap
        self.geneCheckBed = geneCheckBed
        self.annotationBed = annotationBed

    def run(self):
        ################
        # load all relevant information for this comparison
        ################

        #alignment_dict maps alignmentIds (primaryKeys) to pslRow objects
        self.alignmentDict = get_alignment_dict(self.alnPsl)

        #seq_dict is a twoBitFile object for this target genome
        self.seqDict = get_seq_dict(self.seqTwoBit)

        #ref_seq_dict is twoBitFile object for the query genome
        self.refSeqDict = get_seq_dict(self.refSeqTwoBit)

        #transcript attributes maps alignmentIds to gencode basic attributes
        self.attributeDict = get_transcript_attributes(self.gencodeAttributeMap)

        #transcript_dict maps alignmentIds to Transcript objects
        self.transcriptDict = get_transcript_dict(self.geneCheckBed)

        #annotationDict maps transcripts to information in the original annotation
        self.annotationDict = get_transcript_dict(self.annotationBed)

        #alnIds stores all alignment IDs
        self.alnIds = self.alignmentDict.keys()

        ################
        # initialize the sqlite3 database
        ################  

        #find all user-defined classes in the classifiers module
        classifiers = classesInModule(src.classifiers)
        #find all column definitions
        column_definitions = [[x.__name__, x._getType()] for x in classifiers]
        initialize_sql_table(self.db, self.genome, column_definitions, self.primaryKey)

        ################
        # begin classifying
        # write SQL every 5000 records
        ################
        
        vd = {}
        with sql_lib.ExclusiveSqlConnection(self.db) as cur:    
            for aId in self.alnIds:
                vd[aId] = [x().classify(aId, self.alignmentDict, self.refSeqDict, self.seqDict, 
                        self.attributeDict, self.transcriptDict, self.annotationDict) 
                        for x in classifiers]
                if len(vd) % 5000 == 0:
                    insert_row_dict_wrapper(cur, self.genome, self.primaryKey, classifiers, vd)
                    vd = {}


def insert_row_dict_wrapper(cur, genome, primaryKeyColumn, classifiers, value_dict):
    """
    wrapper for sql_lib insertRows where the values are a dict keyed on the primary key.
    classifiers is a list of columns to be filled out.
    """
    column_string = ", ".join([primaryKeyColumn] + [x.__name__ for x in classifiers])

    def value_iter(value_dict):
        for aId, vals in value_dict.iteritems():
            yield [aId] + vals

    sql_lib.insertRows(cur, genome, column_string, len(classifiers) + 1, value_iter(value_dict))


def get_transcript_attributes(gencodeAttributeMap):
    return seq_lib.getTranscriptAttributeDict(gencodeAttributeMap)


def get_transcript_dict(geneCheckBed):
    transcripts = seq_lib.getTranscripts(geneCheckBed)
    return seq_lib.transcriptListToDict(transcripts, noDuplicates=True)


def get_seq_dict(seqTwoBit):
    return seq_lib.readTwoBit(seqTwoBit)


def get_alignment_dict(psl):
    psl_list = psl_lib.readPsl(psl)
    return psl_lib.getPslDict(psl_list, noDuplicates=True)


def initialize_sql_table(db, genome, columns, primaryKey):
    """
    initializes a sql table at <db> based on columns.
    table name will be <genome>, and have a unnique text primary key
    whose column name is <primary_key>.
    columns is a list of pairs [<col_name>,<data_type>] that is each classifer.
    """
    with sql_lib.ExclusiveSqlConnection(db) as cur:
        sql_lib.initializeTable(cur, genome, columns, primaryKey)
