"""
Alignment attributes for comparativeAnnotator.
"""
from sqlalchemy import Float, Integer, String
import comparativeAnnotator.lib.annotation_utils as utils
from pycbio.sys.mathOps import format_ratio

__author__ = "Ian Fiddes"


class AlignmentCoverage(utils.AbstractClassifier):
    """
    Calculates alignment coverage:
    (matches + mismatches + repeat matches) / q_size
    Reports the value as a REAL between 0 and 1
    """
    dtype = Float

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return aln.coverage


class AlignmentIdentity(utils.AbstractClassifier):
    """
    Calculates alignment identity:
    matches / (matches + mismatches + query_insertions)
    Reports the value as a REAL between 0 and 1
    """
    dtype = Float

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return aln.identity


class PercentUnknownBases(utils.AbstractClassifier):
    """
    Calculates the percent of unknown bases in the alignment:
    n_count / q_size
    """
    dtype = Float

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return aln.percent_n


class PercentUnknownCodingBases(utils.AbstractClassifier):
    """
    Calculates the percent of coding bases that are Ns in the transcript
    """
    dtype = Float

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        cds = t.get_cds(tgt_fasta)
        return 100 * format_ratio(cds.count("N"), len(cds))


class NumberIntrons(utils.AbstractClassifier):
    """
    Reports the number of introns for this alignment
    """
    dtype = Integer

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return len(t.intron_intervals)


class NumberMissingOriginalIntrons(utils.AbstractClassifier):
    """
    Does the alignment have all original introns? It can have more (small gaps and such), but it must have all
    original introns. Reports the number of missing introns.
    """
    dtype = Integer

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta, fuzz_distance=5):
        aln_starts_ends = utils.get_adjusted_starts_ends(t, aln)
        count = 0
        for ref_exon in a.exons[1:]:
            r = [aln_start - fuzz_distance <= ref_exon.start <= aln_end + fuzz_distance for aln_start, aln_end in
                 aln_starts_ends]
            if not any(r):
                count += 1
        return count


class TranscriptId(utils.AbstractClassifier):
    """
    Reports the original transcript ID. Used to map between reference tables and target tables.
    """
    dtype = String

    def __call__(self, a, t, aln, ref_aln, ref_fasta, tgt_fasta):
        return a.name