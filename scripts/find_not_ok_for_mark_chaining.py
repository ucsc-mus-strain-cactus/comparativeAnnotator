import argparse
from scripts.plot_functions import *
import lib.sequence_lib as seq_lib
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outDir", required=True, help="output")
    parser.add_argument("--simpleChainDir", required=True, help="directory containing databases")
    parser.add_argument("--allDir", required=True, help="directory containing databases")
    parser.add_argument("--annotationGp", type=str, required=True, help="annotation genePred")
    parser.add_argument("--allGp", type=str, required=True, help="all GP to pull bad transcripts from")
    parser.add_argument("--simpleGp", type=str, required=True, help="simpleChain GP to compare to")
    parser.add_argument("--attributePath", type=str, required=True, help="attribute tsv file")
    return parser.parse_args()


genome = "C57B6NJ"  # the genome we want to find genes that are OK in one set that are not in the other

def parse_details(records):
    for record in records:
        for x in record:
            if x == None:
                continue
            yield x + "\n"


def get_transcript_dict(gp, filter_set):
    transcripts = seq_lib.getGenePredTranscripts(gp)
    d = seq_lib.transcriptListToDict(transcripts, noDuplicates=True)
    r = defaultdict(list)
    for aln_id, rec in d.iteritems():
        tx_id = strip_alignment_numbers(aln_id)
        if tx_id in filter_set:
            r[tx_id].append(rec)
    return r


def get_ids_we_care_about(attr_path, annotation_gp):
    gencode_ids = get_gp_ids(annotation_gp)
    biotype_ids = get_all_ids(attr_path, biotype="protein_coding")
    return gencode_ids & biotype_ids


def get_ok(cur, genome, these_ids):
    return {x for x in transmap_ok(cur, genome, tm_coding_classifiers) if strip_alignment_numbers(x) in these_ids}


def find_corresponding_transcript(a_transcripts, s_transcripts):
    r_map = {}
    for tx_id, a_aln in a_transcripts.iteritems():
        s_aln = s_transcripts[tx_id]
        if len(s_aln) == len(a_aln) == 1:
            r_map[s_aln[0]] = a_aln[0]
        else:
            for s in s_aln:
                for a in a_aln:
                    if s.start == a.start and s.stop == a.stop:
                        assert s not in r_map
                        r_map[s] = a
                        break
    return r_map


def find_not_ok_in_a(s_in_a, s_ok, a_ok):
    to_investigate = []
    for s_tx, a_tx in s_in_a.iteritems():
        if s_tx.name in s_ok and a_tx.name not in a_ok:
            to_investigate.append(a_tx)
    return to_investigate


def write_tx_bed(out_dir, to_investigate):
    with open(os.path.join(out_dir, "not_ok_all_chaining.bed"), "w") as outf:
        outf.write('track name="Transcripts OK in simpleChain and not OK in allChain"\n')
        for t in to_investigate:
            outf.write("\t".join(map(str, t.getBed())) + "\n")


def write_human_readable_classifiers(out_dir, to_investigate, a_con):
    formatted_names = ", ".join(['"' + x.name + '"' for x in to_investigate])
    formatted_classifiers = ", ".join(["AlignmentId"] + tm_coding_classifiers)
    cmd = "SELECT {} FROM C57B6NJ WHERE AlignmentId in ({})".format(formatted_classifiers, formatted_names)
    a_data = pd.read_sql(cmd, a_con)
    failures = defaultdict(list)
    for pos, row in a_data.iterrows():
        for classifier, value in row.iteritems():
            if classifier == "AlignmentId":
                name = value
            elif value == 1:
                failures[name].append(classifier)
    with open(os.path.join(out_dir, "failed_classifiers.tsv"), "w") as outf:
        for name, vals in failures.iteritems():
            vals = ",".join(sorted(vals))
            outf.write("\t".join([name, vals]) + "\n")
    return formatted_names


def write_browser_bed(out_dir, all_dir, formatted_names):
    a_details_con = sql.connect(os.path.join(all_dir, "details.db"))
    a_details_cur = a_details_con.cursor()
    formatted_classifiers = ", ".join(tm_coding_classifiers)
    cmd = "SELECT {} FROM C57B6NJ WHERE AlignmentId in ({})".format(formatted_classifiers, formatted_names)
    recs = a_details_cur.execute(cmd).fetchall()
    with open(os.path.join(out_dir, "failed_classifiers.bed"), "w") as outf:
        outf.write('track name="Classifiers failed in allChain transcripts that were ok in simpleChain"\n')
        for r in parse_details(recs):
            outf.write(r)

def main():
    args = parse_args()
    s_con, s_cur = attach_databases(args.simpleChainDir)
    a_con, a_cur = attach_databases(args.allDir)
    these_ids = get_ids_we_care_about(args.attributePath, args.annotationGp)
    s_ok = get_ok(s_cur, genome, these_ids)
    a_ok = get_ok(a_cur, genome, these_ids)
    a_transcripts = get_transcript_dict(args.allGp, these_ids)
    s_transcripts = get_transcript_dict(args.simpleGp, these_ids)
    s_in_a = find_corresponding_transcript(a_transcripts, s_transcripts)
    to_investigate = find_not_ok_in_a(s_in_a, s_ok, a_ok)
    write_tx_bed(args.outDir, to_investigate)
    formatted_names = write_human_readable_classifiers(args.outDir, to_investigate, a_con)
    write_browser_bed(args.outDir, args.allDir, formatted_names)
    # pull out the BED records for this transcript


if __name__ == "__main__":
    main()