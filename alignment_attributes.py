"""
Alignment attributes for comparativeAnnotator.
"""
from peewee import FloatField, IntegerField
from collections import Counter, defaultdict
from comparativeAnnotator.comp_lib.name_conversions import remove_alignment_number
import comparativeAnnotator.comp_lib.annotation_utils as utils
from pycbio.sys.mathOps import format_ratio

__author__ = "Ian Fiddes"


class AlignmentCoverage(utils.AbstractClassifier):
    """
    Calculates alignment coverage:
    (matches + mismatches + repeat matches) / q_size
    Reports the value as a REAL between 0 and 1
    """
    dtype = FloatField

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return aln.coverage


class AlignmentIdentity(utils.AbstractClassifier):
    """
    Calculates alignment identity:
    matches / (matches + mismatches + query_insertions)
    Reports the value as a REAL between 0 and 1
    """
    dtype = FloatField

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return aln.identity


class PercentUnknownBases(utils.AbstractClassifier):
    """
    Calculates the percent of unknown bases in the alignment:
    n_count / q_size
    """
    dtype = FloatField

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return aln.percent_n


class PercentUnknownCodingBases(utils.AbstractClassifier):
    """
    Calculates the percent of coding bases that are Ns in the transcript
    """
    dtype = FloatField

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        cds = t.get_cds(tgt_fasta)
        return 100 * format_ratio(cds.count("N"), len(cds))


class NumberIntrons(utils.AbstractClassifier):
    """
    Reports the number of introns for this alignment
    """
    dtype = IntegerField

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return len(t.intron_intervals)


class NumberMissingOriginalIntrons(utils.AbstractClassifier):
    """
    Does the alignment have all original introns? It can have more (small gaps and such), but it must have all
    original introns. Reports the number of missing introns.
    """
    dtype = IntegerField

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, fuzz_distance=5):
        aln_starts_ends = utils.get_adjusted_starts_ends(t, aln)
        count = 0
        for ref_exon in a.exons[1:]:
            r = [aln_start - fuzz_distance <= ref_exon.start <= aln_end + fuzz_distance for aln_start, aln_end in
                 aln_starts_ends]
            if not any(r):
                count += 1
        return count


def paralogy(psl_dict):
    """
    This special non-classifier function takes the entire transcript dict and produces counts of paralogy.
    """
    counts = Counter(remove_alignment_number(aln_id) for aln_id in psl_dict.iterkeys())
    results = {}
    for aln_id, psl in psl_dict.iteritems():
        results[aln_id] = counts[remove_alignment_number(aln_id)] - 1
    return results


def highest_cov_aln(psl_dict):
    """
    This special classifier reports whether a given alignment is the highest coverage alignment for that source tx.
    """
    combined_covs = defaultdict(list)
    for aln_id, psl in psl_dict.iteritems():
        tx_id = remove_alignment_number(aln_id)
        combined_covs[tx_id].append([aln_id, psl.coverage])
    best_cov = set()
    for tx_id, vals in combined_covs.iteritems():
        best = sorted(vals, key=lambda (aln_id, cov): -cov)[0][0]
        best_cov.add(best)
    return {aln_id: True if aln_id in best_cov else False for aln_id in psl_dict.iterkeys()}


def highest_ident_aln(psl_dict):
    """
    This special classifier reports whether a given alignment is the highest identity alignment for that source tx.
    """
    combined_idents = defaultdict(list)
    for aln_id, psl in psl_dict.iteritems():
        tx_id = remove_alignment_number(aln_id)
        combined_idents[tx_id].append([aln_id, psl.identity])
    best_ident = set()
    for tx_id, vals in combined_idents.iteritems():
        best = sorted(vals, key=lambda (aln_id, ident): -ident)[0][0]
        best_ident.add(best)
    return {aln_id: True if aln_id in best_ident else False for aln_id in psl_dict.iterkeys()}


def best_overall_aln(psl_dict, cov_weight=0.25, ident_weight=0.75):
    """
    This special classifier reports whether a given alignment has the highest sum of coverage + ident, weighted
    """
    assert cov_weight + ident_weight == 1.0
    combined_stats = defaultdict(list)
    for aln_id, psl in psl_dict.iteritems():
        tx_id = remove_alignment_number(aln_id)
        combined_stats[tx_id].append([aln_id, ident_weight * psl.identity + cov_weight * psl.coverage])
    best_ident = set()
    for tx_id, vals in combined_stats.iteritems():
        best = sorted(vals, key=lambda (aln_id, s): -s)[0][0]
        best_ident.add(best)
    return {aln_id: True if aln_id in best_ident else False for aln_id in psl_dict.iterkeys()}
