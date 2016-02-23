from collections import defaultdict
from name_conversions import strip_alignment_numbers
from pycbio.sys.sqliteOps import get_query_ids, get_query_dict, get_non_unique_query_dict


def get_biotype_aln_ids(cur, genome, biotype, filter_chroms=None):
    """
    Returns a set of aln_ids which are members of the biotype.
    """
    query = "SELECT AlignmentId FROM attributes.'{0}' WHERE TranscriptType = '{1}' AND GeneType = '{1}'"
    if filter_chroms is not None:
        for filter_chrom in filter_chroms:
            query += " AND refChrom != '{}'".format(filter_chrom)
    query = query.format(genome, biotype)
    return get_query_ids(cur, query)


def get_biotype_ids(cur, ref_genome, biotype, mode="Transcript", filter_chroms=None):
    """
    Returns a set of transcript_ids which are members of the biotype. Category should be Transcript or Gene.
    In gene mode returns the gene_id instead of transcript_id.
    """
    assert mode in ["Transcript", "Gene"]
    query = ("SELECT {0}Id FROM main.'{1}' JOIN attributes.'{1}' USING (TranscriptId) WHERE "
             "TranscriptType = '{2}' AND GeneType = '{2}'")
    if filter_chroms is not None:
        for filter_chrom in filter_chroms:
            query += " AND refChrom != '{}'".format(filter_chrom)
    query = query.format(mode, ref_genome, biotype)
    return get_query_ids(cur, query)


def get_all_biotypes(cur, ref_genome, gene_level=False):
    """
    Returns all biotypes in the attribute database.
    """
    r = "Gene" if gene_level else "Transcript"
    query = "SELECT {}Type FROM attributes.'{}'".format(r, ref_genome)
    return get_query_ids(cur, query)


def get_ids_by_chromosome(cur, genome, chromosomes=("Y", "chrY")):
    """
    Returns all AlignmentIds in the attribute database whose source chromosome is in chromomsomes
    """
    base_query = "SELECT AlignmentId FROM attributes.'{}' WHERE {}"
    filt = ["SourceChrom = '{}'".format(x) for x in chromosomes]
    filt = " OR ".join(filt)
    query = base_query.format(genome, filt)
    return get_query_ids(cur, query)


def get_transcript_gene_map(cur, ref_genome, biotype=None, filter_chroms=None):
    """
    Returns a dictionary mapping transcript IDs to gene IDs
    """
    if biotype is None:
        query = "SELECT transcriptId,geneId FROM attributes.'{0}'"
    else:
        query = "SELECT transcriptId,geneId FROM attributes.'{0}' WHERE transcriptType = '{1}' AND geneType = '{1}'"
    if filter_chroms is not None:
        add = "AND".join([" refChrom != '{}' ".format(filter_chrom) for filter_chrom in filter_chroms])
        if biotype is None:
            query += " WHERE "
        else:
            query += " AND "
        query += add
    query = query.format(ref_genome, biotype)
    return get_query_dict(cur, query)


def get_gene_transcript_map(cur, ref_genome, biotype=None, filter_chroms=None):
    """
    Returns a dictionary mapping all gene IDs to their respective transcripts
    """
    if biotype is None:
        query = "SELECT geneId,transcriptId FROM attributes.'{0}'"
    else:
        query = "SELECT geneId,transcriptId FROM attributes.'{0}' WHERE transcriptType = '{1}' AND geneType = '{1}'"
    if filter_chroms is not None:
        add = "AND".join([" refChrom != '{}' ".format(filter_chrom) for filter_chrom in filter_chroms])
        if biotype is None:
            query += " WHERE "
        else:
            query += " AND "
        query += add
    query = query.format(ref_genome, biotype)
    return get_non_unique_query_dict(cur, query)


def get_gene_biotype_map(cur, ref_genome):
    """
    Returns a dictionary mapping all gene IDs to their respective biotypes
    """
    query = "SELECT geneId,geneType FROM attributes.'{}'".format(ref_genome)
    return get_query_dict(cur, query)


def get_transcript_biotype_map(cur, ref_genome):
    """
    Returns a dictionary mapping all transcript IDs to their respective biotypes
    """
    query = "SELECT transcriptId,transcriptType FROM attributes.'{}'".format(ref_genome)
    return get_query_dict(cur, query)


def get_fail_passing_excel_ids(cur, ref_genome, genome, biotype, best_cov_only=True, filter_chroms=None,
                           highest_cov_dict=None):
    """
    Returns the IDs categorized as fail, passing_specific, excellent. You can set the best_cov_only flag to only report
    those transcripts with highest coverage. You can also excellent a premade highest_cov_dict to save computation time.
    """
    passing_query = etc.config.transMapEval(ref_genome, genome, biotype, passing=True)
    excellent_query = etc.config.transMapEval(ref_genome, genome, biotype, passing=False)
    if best_cov_only is True:
        if highest_cov_dict is None:
            best_covs = highest_cov_aln(cur, genome)
        else:
            best_covs = highest_cov_dict[genome]
        best_ids = set(zip(*best_covs.itervalues())[0])
    else:
        class Universe:
            def __contains__(_, __):
                return True
        best_ids = Universe()
    passing_ids = {x for x in get_query_ids(cur, passing_query) if x in best_ids}
    excellent_ids = {x for x in get_query_ids(cur, excellent_query) if x in best_ids}
    all_ids = {x for x in get_biotype_aln_ids(cur, genome, biotype) if x in best_ids}
    fail_ids = all_ids - (passing_ids | excellent_ids)
    passing_specific_ids = passing_ids - excellent_ids
    return fail_ids, passing_specific_ids, excellent_ids


def collapse_details_dict(details_dict):
    """
    Collapses a details dict into a string so that the database record is a string. Properly accounts for the three
    possibilities: the record is a single BED entry (a string), the record is a list of BED entries, or the record
    is empty but has been accessed (because the dict is a defaultdict)
    """
    collapsed = {}
    for aln_id, rec in details_dict.iteritems():
        if len(rec) == 0:  # to deal with the downsides of defaultdict
            continue
        if isinstance(rec[0], list):  # hack to determine if this is a list of lists
            rec = "\n".join(["\t".join(map(str, x)) for x in rec])
        else:
            rec = "\t".join(map(str, rec))
        collapsed[aln_id] = "".join([rec, "\n"])
    return collapsed


def run_transmap_eval(cur, genome, biotype, trans_map_eval, passing=True):
    """
    Convenience wrapper for getting all Pass/Excellent transcripts for transMap
    """
    query = trans_map_eval(genome, biotype, passing)
    return get_query_ids(cur, query)


def get_stats(cur, genome, mode="transMap", filter_chroms=None):
    """
    Returns a dictionary mapping each aln_id to [aln_id, %ID, %COV].
    We replace all NULL values with zero to make things easier in consensus finding.
    If filter_chroms is set, will filter out transcripts that come from those transcripts. However, at this point,
    if this is done in augustus mode nothing will happen.
    """
    if mode == "transMap":
        query = ("SELECT AlignmentId,IFNULL(AlignmentCoverage, 0),IFNULL(AlignmentIdentity, 0) "
                 "FROM attributes.'{}'")
        if filter_chroms is not None:
            add = "AND".join([" sourceChrom != '{}' ".format(filter_chrom) for filter_chrom in filter_chroms])
            query += " WHERE "
            query += add
    else:
        query = ("SELECT AugustusAlignmentId,IFNULL(AlignmentCoverage, 0),IFNULL(AlignmentIdentity, 0) "
                 "FROM augustus_attributes.'{}'")
    return get_query_dict(cur, query.format(genome))


def highest_cov_aln(cur, genome, filter_chroms=None):
    """
    Returns the set of alignment IDs that represent the best alignment for each source transcript (that mapped over)
    Best is defined as highest %COV. Also reports the associated coverage and identity values.
    """
    tm_stats = get_stats(cur, genome, mode="transMap", filter_chroms=filter_chroms)
    combined_covs = defaultdict(list)
    for aln_id, (cov, ident) in tm_stats.iteritems():
        tx_id = strip_alignment_numbers(aln_id)
        combined_covs[tx_id].append([aln_id, cov, ident])
    best_cov = {}
    for tx_id, vals in combined_covs.iteritems():
        best_cov[tx_id] = sorted(vals, key=lambda (aln_id, cov, ident): (-cov, -ident))[0]
    return best_cov


def get_highest_cov_alns(cur, genomes, filter_chroms=None):
    """
    Dictionary mapping each genome to a dictionary reporting each highest coverage alignment and its metrics
    """
    return {genome: highest_cov_aln(cur, genome, filter_chroms=filter_chroms) for genome in genomes}
