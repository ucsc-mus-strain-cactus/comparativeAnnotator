import sys
import os
import re
import argparse
import itertools
from collections import defaultdict
import lib.sequence_lib as seq_lib
import lib.sqlite_lib as sql_lib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refGp", required=True)
    parser.add_argument("--refFasta", required=True)
    parser.add_argument('--refGenome', type=str, required=True)
    parser.add_argument('--gencodeAttributeMap', required=True)
    parser.add_argument('--primaryKeyColumn', type=str, default="ensId")
    parser.add_argument('--outDir', type=str, required=True)
    return parser.parse_args()

rgb = "152,156,45"


def get_transcript_dict(gp_path):
    """
    Loads the reference genePred as Transcript objects
    """
    transcripts = seq_lib.getGenePredTranscripts(gp_path)
    transcript_dict = seq_lib.transcriptListToDict(transcripts, noDuplicates=True)
    return transcript_dict


def transcript_iterator(transcript_dict):
    """
    Convenience function for iterating over a dictionary of Transcript objects
    """
    for ens_id, t in transcript_dict.iteritems():
        yield ens_id, t


def AbutsUnknownBases(transcript_dict, seq_dict, short_intron_size=30, distance=10):
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        intervals = [[t.exonIntervals[0].start - distance, t.exonIntervals[0].start]]
        for intron in t.intronIntervals:
            if len(intron) > short_intron_size:
                intervals.append([intron.start, intron.start + distance])
        intervals.append([t.exonIntervals[-1].stop, t.exonIntervals[-1].stop + distance])
        for start, stop in intervals:
            seq = seq_dict[t.chromosome][start:stop]
            if "N" in seq:
                classify_dict[ens_id] = 1
                details_dict[ens_id].append(t.getBed(rgb, sys._getframe().f_code.co_name))
                break
        if ens_id not in classify_dict:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def StartOutOfFrame(transcript_dict):
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        if t.getCdsLength() <= 75:
            continue
        t_frames = [x for x in t.exonFrames if x != -1]
        if t.strand is True and t_frames[0] != 0 or t.strand is False and t_frames[-1] != 0:
            classify_dict[ens_id] = 1
            details_dict[ens_id] = seq_lib.cdsCoordinateToBed(t, 0, 3, rgb, sys._getframe().f_code.co_name)
            continue
        classify_dict[ens_id] = 0
    return classify_dict, details_dict


def BadFrame(transcript_dict):
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        if t.getCdsLength() <= 75:
            continue
        if t.getCdsLength() % 3 != 0:
                classify_dict[ens_id] = 1
                details_dict[ens_id].append(t.getBed(rgb, sys._getframe().f_code.co_name))
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def BeginStart(transcript_dict, seq_dict):
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        seq = t.getCds(seq_dict)
        s = len(s)
        if s <= 75:
            continue
        if seq[:3] != "ATG":
            classify_dict[ens_id] = 1
            details_dict[ens_id] = seq_lib.cdsCoordinateToBed(t, 0, 3, rgb, sys._getframe().f_code.co_name)
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def EndStop(transcript_dict, seq_dict):
    stop_codons = ('TAA', 'TGA', 'TAG')
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        seq = t.getCds(seq_dict)
        s = len(s)
        if s <= 75:
            continue
        if seq[-3:] not in stop_codons:
            classify_dict[ens_id] = 1
            details_dict[ens_id] = seq_lib.cdsCoordinateToBed(t, s - 3, s, rgb,  sys._getframe().f_code.co_name)
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def base_gap(transcript_dict, seq_dict, inequality, short_intron_size, skip_n):
    classify_dict = {}
    details_dict = defaultdict(list)
    for ens_id, t in transcript_iterator(transcript_dict):
        for intron in t.intronIntervals:
            if len(intron) >= short_intron_size:
                continue
            elif skip_n and "N" in intron.getSequence(seq_dict):
                continue
            elif inequality(intron, t):
                continue
            classify_dict[ens_id] = 1
            details_dict[ens_id].append(seq_lib.intervalToBed(t, intron, rgb, sys._getframe().f_code.co_name))
        if ens_id not in classify_dict:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def is_cds(intron, t):
    return not (intron.start >= t.thickStart and intron.stop < t.thickStop)


def is_not_cds(intron, t):
    return intron.start >= t.thickStart and intron.stop < t.thickStop


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
            seq = intron.getSequence(seq_dict, strand=True)
            donor, acceptor = seq[:2], seq[-2:]
            if donor not in splice_sites or splice_sites[donor] != acceptor:
                classify_dict[ens_id] = 1
                details_dict[ens_id].append(seq_lib.spliceIntronIntervalToBed(t, intron, rgb,
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
        cds = t.getCds(seq_dict)
        offset = seq_lib.findOffset(t.exonFrames, t.strand)
        for i, codon in seq_lib.readCodonsWithPosition(cds, offset, skip_last=True):
            amino_acid = seq_lib.codonToAminoAcid(codon)
            if amino_acid == "*":
                classify_dict[ens_id] = 1
                details_dict[ens_id].append(seq_lib.cdsCoordinateToBed(t, i, i + 3, rgb,
                                                                       sys._getframe().f_code.co_name))
        if ens_id not in classify_dict:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def ShortCds(transcript_dict, cds_cutoff=75):
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        if t.getCdsLength() < 3:
            continue
        elif t.getCdsLength() <= cds_cutoff:
            classify_dict[ens_id] = 1
            details_dict[ens_id] = seq_lib.transcriptToBed(t, rgb, sys._getframe().f_code.co_name)
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def unknown_base(transcript_dict, seq_dict, r, cds):
    classify_dict = {}
    details_dict = {}
    for ens_id, t in transcript_iterator(transcript_dict):
        if cds is True:
            s = t.getCds(seq_dict)
            tmp = [seq_lib.cdsCoordinateToBed(t, m.start(), m.end(), rgb, sys._getframe().f_code.co_name) for m in
                   re.finditer(r, s)]
        else:
            s = t.getMRna(seq_dict)
            tmp = [seq_lib.transcriptCoordinateToBed(t, m.start(), m.end(), rgb, sys._getframe().f_code.co_name)
                   for m in re.finditer(r, s)]
        if len(tmp) > 0:
            details_dict[ens_id] = tmp
            classify_dict[ens_id] = 1
        else:
            classify_dict[ens_id] = 0
    return classify_dict, details_dict


def UnknownBases(transcript_dict, seq_dict):
    r = re.compile("[atgcATGC][N]{1,10}[atgcATGC]")
    return unknown_base(transcript_dict, seq_dict, r, cds=False)


def UnknownCdsBases(transcript_dict, seq_dict):
    r = re.compile("[atgcATGC][N]{1,10}[atgcATGC]")
    return unknown_base(transcript_dict, seq_dict, r, cds=True)


def ScaffoldGap(transcript_dict, seq_dict):
    r = re.compile("[atgcATGC][N]{11,}[atgcATGC]")
    return unknown_base(transcript_dict, seq_dict, r, cds=False)


def invert_dict(d):
    for a, b in d.iteritems():
        yield b, a


def get_ids(gp):
    return {x.split()[0] for x in open(gp)}


def initialize_db(db, classifiers, ref_genome, primary_key, ref_gp, data_type="TEXT"):
    ens_ids = get_ids(ref_gp)
    assert data_type in ["TEXT", "INTEGER"]
    with sql_lib.ExclusiveSqlConnection(db) as cur:
        sql_lib.initializeTable(cur, ref_genome, classifiers, primary_key)
        sql_lib.insertRows(cur, ref_genome, primary_key, [primary_key], itertools.izip_longest(ens_ids, [None]))


def main():
    args = parse_args()
    fns = [AbutsUnknownBases, StartOutOfFrame, BadFrame, BeginStart, EndStop, CdsGap, UtrGap, UnknownGap,
           CdsNonCanonSplice, CdsUnknownSplice, UtrNonCanonSplice, UtrUnknownSplice, InFrameStop, ShortCds,
           UnknownBases, UnknownCdsBases, ScaffoldGap]
    classify_dicts = {}
    details_dicts = {}
    fn_args = {"transcript_dict": get_transcript_dict(args.refGp), "seq_dict": seq_lib.getSequenceDict(args.refFasta)}
    for fn in fns:
        classify_dicts[fn.__name__], details_dicts[fn.__name__] = fn(**fn_args)
    details_db = os.path.join(args.outDir, "details.db")
    classify_db = os.path.join(args.outDir, "classify.db")
    classifier_names = [x.__name__ for x in fns]
    initialize_db(details_db, classifier_names, args.refGenome, args.primaryKeyColumn, args.refGp, data_type="TEXT")
    initialize_db(classify_db, classifier_names, args.refGenome, args.primaryKeyColumn, args.refGp, data_type="INTEGER")
    with sql_lib.ExclusiveSqlConnection(classify_db) as cur:
        for column, value_dict in classify_dicts.iteritems():
            sql_lib.updateRows(cur, args.refGenome, args.primaryKeyColumn, column, invert_dict(value_dict))
    with sql_lib.ExclusiveSqlConnection(details_db) as cur:
        for column, value_dict in details_dicts.iteritems():
            sql_lib.updateRows(cur, args.refGenome, args.primaryKeyColumn, column, invert_dict(value_dict))
