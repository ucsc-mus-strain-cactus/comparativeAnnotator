"""
This script runs a specific subset of classifiers on the reference. This script is designed to run without jobTree
and can be applied to any genePred (is not comparative).
"""

import sys
import re
import argparse
import itertools
from collections import defaultdict
import lib.seq_lib as seq_lib
import lib.sql_lib as sql_lib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refGp", required=True)
    parser.add_argument("--refFasta", required=True)
    parser.add_argument('--refGenome', type=str, required=True)
    parser.add_argument('--outDir', type=str, required=True)
    return parser.parse_args()

rgb = "152,156,45"


def get_transcript_dict(gp_path):
    """
    Loads the reference genePred as Transcript objects
    """
    transcripts = seq_lib.get_gene_pred_transcripts(gp_path)
    transcript_dict = seq_lib.transcript_list_to_dict(transcripts)
    return transcript_dict


def transcript_iterator(transcript_dict):
    """
    Convenience function for iterating over a dictionary of Transcript objects
    """
    for ens_id, t in transcript_dict.iteritems():
        yield ens_id, t


def AbutsUnknownBases(transcript_dict, seq_dict, short_intron_size=30, distance=2):
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        intervals = [[t.exon_intervals[0].start - distance, t.exon_intervals[0].start]]
        for intron in t.intron_intervals:
            if len(intron) > short_intron_size:
                intervals.append([intron.start, intron.start + distance])
        intervals.append([t.exon_intervals[-1].stop, t.exon_intervals[-1].stop + distance])
        for start, stop in intervals:
            seq = seq_dict[t.chromosome][start:stop]
            if "N" in seq:
                classify_dict[ens_id] = 1
                details_dict[ens_id].append(t.get_bed(rgb, sys._getframe().f_code.co_name))
                break
        if ens_id not in classify_dict:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def StartOutOfFrame(transcript_dict, seq_dict):
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        if t.getCdsLength() <= 75:
            continue
        t_frames = [x for x in t.exon_frames if x != -1]
        if t.strand is True and t_frames[0] != 0 or t.strand is False and t_frames[-1] != 0:
            classify_dict[ens_id] = 1
            details_dict[ens_id] = seq_lib.cds_coordinate_to_bed(t, 0, 3, rgb, sys._getframe().f_code.co_name)
            continue
        classify_dict[ens_id] = 0
    return classify_dict, details_dict


def BadFrame(transcript_dict, seq_dict):
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        if t.get_cds_length() <= 75:
            continue
        if t.get_cds_length() % 3 != 0:
                classify_dict[ens_id] = 1
                details_dict[ens_id].append(t.get_bed(rgb, sys._getframe().f_code.co_name))
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def BeginStart(transcript_dict, seq_dict):
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        seq = t.get_cds(seq_dict)
        if len(seq) <= 75:
            continue
        if seq[:3] != "ATG":
            classify_dict[ens_id] = 1
            details_dict[ens_id] = seq_lib.cds_coordinate_to_bed(t, 0, 3, rgb, sys._getframe().f_code.co_name)
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def EndStop(transcript_dict, seq_dict):
    stop_codons = ('TAA', 'TGA', 'TAG')
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        seq = t.get_cds(seq_dict)
        s = len(seq)
        if s <= 75:
            continue
        if seq[-3:] not in stop_codons:
            classify_dict[ens_id] = 1
            details_dict[ens_id] = seq_lib.cds_coordinate_to_bed(t, s - 3, s, rgb,  sys._getframe().f_code.co_name)
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def base_gap(transcript_dict, seq_dict, inequality, short_intron_size, skip_n):
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        for intron in t.intron_intervals:
            if len(intron) >= short_intron_size:
                continue
            elif skip_n and "N" in intron.get_sequence(seq_dict):
                continue
            elif inequality(intron, t):
                continue
            classify_dict[ens_id] = 1
            details_dict[ens_id].append(seq_lib.interval_to_bed(t, intron, rgb, sys._getframe().f_code.co_name))
        if ens_id not in classify_dict:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def is_cds(intron, t):
    return not (intron.start >= t.thick_start and intron.stop < t.thick_stop)


def is_not_cds(intron, t):
    return intron.start >= t.thick_start and intron.stop < t.thick_stop


def dummy_inequality(intron, t):
    return True


def CdsGap(transcript_dict, seq_dict, short_intron_size=30):
    return base_gap(transcript_dict, seq_dict, is_cds, short_intron_size, skip_n=True)


def UtrGap(transcript_dict, seq_dict, short_intron_size=30):
    return base_gap(transcript_dict, seq_dict, is_not_cds, short_intron_size, skip_n=True)


def UnknownGap(transcript_dict, seq_dict, short_intron_size=30):
    return base_gap(transcript_dict, seq_dict, dummy_inequality, short_intron_size, skip_n=False)


def base_splice(transcript_dict, seq_dict, inequality, short_intron_size, splice_sites):
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        for intron in t.intronIntervals:
            if len(intron) <= short_intron_size:
                continue
            elif inequality(intron, t):
                continue
            seq = intron.get_sequence(seq_dict, strand=True)
            donor, acceptor = seq[:2], seq[-2:]
            if donor not in splice_sites or splice_sites[donor] != acceptor:
                classify_dict[ens_id] = 1
                details_dict[ens_id].append(seq_lib.splice_intron_interval_to_bed(t, intron, rgb,
                                                                              sys._getframe().f_code.co_name))
        if ens_id not in classify_dict:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def CdsNonCanonSplice(transcript_dict, seq_dict, short_intron_size=30):
    canonical = {"GT": "AG"}
    return base_splice(transcript_dict, seq_dict, is_cds, short_intron_size, canonical)


def CdsUnknownSplice(transcript_dict, seq_dict, short_intron_size=30):
    non_canonical = {"GT": "AG", "GC": "AG", "AT": "AC"}
    return base_splice(transcript_dict, seq_dict, is_cds, short_intron_size, non_canonical)


def UtrNonCanonSplice(transcript_dict, seq_dict, short_intron_size=30):
    canonical = {"GT": "AG"}
    return base_splice(transcript_dict, seq_dict, is_not_cds, short_intron_size, canonical)


def UtrUnknownSplice(transcript_dict, seq_dict, short_intron_size=30):
    non_canonical = {"GT": "AG", "GC": "AG", "AT": "AC"}
    return base_splice(transcript_dict, seq_dict, is_not_cds, short_intron_size, non_canonical)


def InFrameStop(transcript_dict, seq_dict):
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        cds = t.get_cds(seq_dict)
        offset = seq_lib.find_offset(t.exonFrames, t.strand)
        for i, codon in seq_lib.read_codons_with_position(cds, offset, skip_last=True):
            amino_acid = seq_lib.codon_to_amino_acid(codon)
            if amino_acid == "*":
                classify_dict[ens_id] = 1
                details_dict[ens_id].append(seq_lib.cds_coordinate_to_bed(t, i, i + 3, rgb,
                                                                       sys._getframe().f_code.co_name))
        if ens_id not in classify_dict:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def ShortCds(transcript_dict, seq_dict, cds_cutoff=75):
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        if t.getCdsLength() < 3:
            continue
        elif t.getCdsLength() <= cds_cutoff:
            classify_dict[ens_id] = 1
            details_dict[ens_id] = seq_lib.transcript_to_bed(t, rgb, sys._getframe().f_code.co_name)
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def unknown_base(transcript_dict, seq_dict, r, cds):
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        if cds is True:
            s = t.get_cds(seq_dict)
            tmp = [seq_lib.cds_coordinate_to_bed(t, m.start(), m.end(), rgb, sys._getframe().f_code.co_name) for m in
                   re.finditer(r, s)]
        else:
            s = t.get_mrna(seq_dict)
            tmp = [seq_lib.transcript_coordinate_to_bed(t, m.start(), m.end(), rgb, sys._getframe().f_code.co_name)
                   for m in re.finditer(r, s)]
        if len(tmp) > 0:
            details_dict[ens_id] = tmp
            classify_dict[ens_id] = 1
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def UnknownBases(transcript_dict, seq_dict):
    r = re.compile("[atgcATGC][N]+[atgcATGC]")
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        s = t.get_mrna(seq_dict)
        tmp = [seq_lib.transcript_coordinate_to_bed(t, m.start() + 1, m.end() - 1, rgb, sys._getframe().f_code.co_name)
               for m in re.finditer(r, s)]
        if len(tmp) > 0:
            details_dict[ens_id] = tmp
            classify_dict[ens_id] = 1
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def invert_dict(d):
    for a, b in d.iteritems():
        yield b, a


def get_ids(gp):
    return {x.split()[0] for x in open(gp)}


def main():
    args = parse_args()
    classifiers = [AbutsUnknownBases, StartOutOfFrame, BadFrame, BeginStart, EndStop, CdsGap, UtrGap, UnknownGap,
                   CdsNonCanonSplice, CdsUnknownSplice, UtrNonCanonSplice, UtrUnknownSplice, InFrameStop, ShortCds,
                   UnknownBases]
    classify_dicts = {}
    details_dicts = {}
    fn_args = {"transcript_dict": get_transcript_dict(args.refGp), "seq_dict": seq_lib.get_sequence_dict(args.refFasta)}
    for fn in classifiers:
        classify_dicts[fn.__name__], details_dicts[fn.__name__] = fn(**fn_args)
    for data_dict, database_path in itertools.izip(*[[classify_dicts, details_dicts], ["classify.db", "details.db"]]):
        sql_lib.write_dict(data_dict, database_path, args.refGenome)


if __name__ == "__main__":
    main()