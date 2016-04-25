"""
jobTree wrapper for AugustusTMR.
TODO: remove sonLib, use pycbio
TODO: re-write to analyze chunks of transcripts, and merge in each of those jobs.
"""

import os
import shutil
import sqlite3 as sql
from pyfasta import Fasta
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import system, popenCatch, getRandomAlphaNumericString, catFiles, TempFileTree
from pycbio.bio.transcripts import GenePredTranscript
from pycbio.sys.fileOps import tmpFileGet, ensureDir
from pycbio.sys.dataOps import grouper


#####
# Below are hard coded commands and SQL queries to run Augustus
#####

hints_db_query = '''select source,typename,start,end,score,strand,frame,priority,grp,mult,esource
      FROM hints,featuretypes
      WHERE hints.speciesid IN (SELECT speciesid FROM speciesnames WHERE speciesname="{genome}")
             AND seqnr IN (SELECT seqnr FROM seqnames WHERE speciesid IN (SELECT speciesid FROM speciesnames WHERE speciesname="{genome}") AND seqname="{chrom}")
             AND start >= {start} AND end <= {stop} AND typeid=type'''

augustus_cmd = ("{augustus} {fasta} --predictionStart=-{start} --predictionEnd=-{start} --extrinsicCfgFile={cfg} "
                "--hintsfile={hints} --UTR=on --alternatives-from-evidence=0 --species=human "
                "--allow_hinted_splicesites=atac --protein=0 --/augustus/verbosity=1 --softmasking=1 "
                "--outfile=/dev/stdout")


def attach_database(path):
    con = sql.connect(path)
    cur = con.cursor()
    return con, cur


def get_transmap_hints(gp_string, tm_2_hints_cmd):
    """
    Uses transMap2hints.pl to create a hints file from a genePred
    """
    return popenCatch(tm_2_hints_cmd, stdinString=gp_string)


def parse_rnaseq_query(query, chromosome, priority=3):
    """
    Parses the results from the sqlite database into a gff-like format
    """
    for source, typename, start, end, score, strand, frame, _, grp, mult, esource in query:
        tags = "pri={};src={};mult={};".format(priority, esource, mult)
        yield "\t".join(map(str, [chromosome, source, typename, start, end, score, ".", ".", tags])) + "\n"


def get_rnaseq_hints(genome, chrom, start, stop, hints_db):
    """
    Extracts the RNAseq hints from the database
    """
    this_db_query = hints_db_query.format(genome=genome, chrom=chrom, start=start, stop=stop)
    con, cur = attach_database(hints_db)
    query = cur.execute(this_db_query)
    rnaseq_hint = "".join(list(parse_rnaseq_query(query, chrom)))
    return rnaseq_hint


def write_hint_fasta(hint, seq, chrom, tmp_dir):
    """
    Writes the hints and the seq to a file to be used by Augustus.
    """
    hint_f = os.path.join(tmp_dir, getRandomAlphaNumericString(10) + ".gff")
    seq_f = os.path.join(tmp_dir, getRandomAlphaNumericString(10) + ".fa")
    with open(hint_f, "w") as hint_fh, open(seq_f, "w") as seq_fh:
        hint_fh.write(hint)
        seq_fh.write(">{}\n{}\n".format(chrom, seq))
    return hint_f, seq_f


def rename_transcripts(transcripts, cfg_version, name):
    """
    Renames overlapping transcripts augIX-ID, where X is the index of the extrinsic.cfg file, e.g. 1 or 2 and where
    ID is the transmap alignment ID, use augIX-n-ID if more than 1 transcript overlaps the alignment
    """
    name_map = {}
    for i, x in enumerate(transcripts):
        if i > 0:
            name_map[x] = "augI{}-{}-{}".format(cfg_version, i + 1, name)
        else:
            name_map[x] = "augI{}-{}".format(cfg_version, name)
    return name_map


def write_augustus(r, name_map, outf_h):
    """
    Writes the results of AugustusTMR to a file.
    """
    features = {"exon", "CDS", "start_codon", "stop_codon", "tts", "tss"}
    for x in r:
        if x.startswith("#"):
            continue
        if "AUGUSTUS" in x:
            x = x.split("\t")
            if x[2] in features:
                t = x[-1].split()
                n = t[-3].split('"')[1]
                if n not in name_map:
                    continue  # skip transcripts filtered out previously
                t[-1] = t[-3] = '"{}";'.format(name_map[n])
                t = " ".join(t)
                x[-1] = t
                outf_h.write("\t".join(map(str, x)) + "\n")


def run_augustus(hint_f, seq_f, name, start, cfg_version, cfg_path, outf_h, gp, augustus_bin):
    """
    Runs Augustus for each cfg/gp_string pair.
    """
    cmd = augustus_cmd.format(fasta=seq_f, start=start, cfg=cfg_path, hints=hint_f, augustus=augustus_bin)
    r = popenCatch(cmd)
    r = r.split("\n")
    # extract only the transcript lines
    l = [x.split() for x in r if "\ttranscript\t" in x]
    # filter out transcripts that do not overlap the alignment range
    transcripts = [x[-1] for x in l if not (int(x[4]) < gp.start or int(x[3]) > gp.stop)]
    # if we lose everything, stop here
    if len(transcripts) > 0:
        # rename transcript based on cfg version, and make names unique
        name_map = rename_transcripts(transcripts, cfg_version, name)
        # write this to a shared location where we will combine later
        write_augustus(r, name_map, outf_h)


def transmap_2_aug(target, gp_strings, args, outf, augustus_bin):
    """
    Runs Augustus on one individual genePred string. Augustus is ran with each cfg file in cfgs
    """
    fasta = Fasta(args.fasta)
    chrom_sizes = {x.split()[0]: x.split()[1] for x in open(args.chrom_sizes)}
    with open(outf, 'w') as outf_h:
        for gp_string in gp_strings:
            gp = GenePredTranscript(gp_string.rstrip().split("\t"))
            # ignore genes with no coding region or longer than max_gene_size
            if not (gp.thick_start >= gp.thick_stop or gp.stop - gp.start > args.max_gene_size):
                chrom = gp.chromosome
                start = max(gp.start - args.padding, 0)
                stop = min(gp.stop + args.padding, chrom_sizes[chrom])
                tm_hint = get_transmap_hints(gp_string, args.tm_2_hints_cmd)
                if args.hints_db is not None:
                    rnaseq_hint = get_rnaseq_hints(args.genome, chrom, start, stop, args.hints_db)
                    hint = "".join([tm_hint, rnaseq_hint])
                else:
                    hint = tm_hint
                seq = fasta[chrom][start:stop]
                hint_f, seq_f = write_hint_fasta(hint, seq, chrom, target.getGlobalTempDir())
                for cfg_version, cfg_path in args.cfgs.iteritems():
                    run_augustus(hint_f, seq_f, gp.name, start, cfg_version, cfg_path, outf_h, gp, augustus_bin)
                # delete the seq and hint file. This makes the final tree cleanup not take so long.
                os.remove(hint_f)
                os.remove(seq_f)


def cat(target, output_gtf, unsorted_tmp_file, out_file_tree):
    """
    Concatenates all of the results into one big GTF, and sorts it by chromosome/pos
    """
    catFiles(out_file_tree.listFiles(), unsorted_tmp_file)
    tmp = tmpFileGet()
    system("sort -k1,1 -k4,4n {} > {}".format(unsorted_tmp_file, tmp))
    ensureDir(os.path.dirname(output_gtf))
    try:
        os.rename(tmp, output_gtf)
    except OSError:
        shutil.copy(tmp, output_gtf)
        os.remove(tmp)


def augustus_tmr_wrapper(target, args):
    """
    Produces one jobTree target per genePred entry.
    """
    # create a file tree in the global output directory. This tree will store the gtf created by each Augustus instance
    out_file_tree = TempFileTree(target.getGlobalTempDir())
    # this file will be where we reduce the final results to before sorting
    unsorted_tmp_file = os.path.join(target.getGlobalTempDir(), getRandomAlphaNumericString(10))
    for line in grouper(open(args.input_gp), 10):
        outf = out_file_tree.getTempFile()
        target.addChildTargetFn(transmap_2_aug, memory=8 * (1024 ** 3),
                                args=[line, args, outf, args.augustus_bin])
    target.setFollowOnTargetFn(cat, args=[args.out_gtf, unsorted_tmp_file, out_file_tree])


def augustus_tmr(args):
    i = Stack(Target.makeTargetFn(augustus_tmr_wrapper, args=[args])).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
