"""
This program takes a series of BAM files and converts them to a Augustus hints database.
"""
import sys
import os
import pysam
import pyfasta
import time
import argparse
import itertools
import subprocess
from collections import defaultdict
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from pycbio.sys.mathOps import format_ratio
from pycbio.sys.fileOps import tmpFileGet, ensureDir, ensureFileDir, atomicInstall
from pycbio.sys.procOps import runProc, runProcCode
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from sonLib.bioio import getRandomAlphaNumericString, catFiles, TempFileTree


def get_tmp(target, global_dir=False, name=None):
    """
    Wrapper functionality for getting temp dirs from jobTree.
    """
    prefix = getRandomAlphaNumericString(10)
    name = "".join([prefix, name]) if name is not None else prefix
    if global_dir is False:
        return os.path.join(target.getLocalTempDir(), name)
    else:
        return os.path.join(target.getGlobalTempDir(), name)


def bam_is_paired(path, num_reads=100000, paired_cutoff=0.75):
    """
    Infers the paired-ness of a bam file.
    """
    sam = pysam.Samfile(path)
    count = 0
    for i, rec in enumerate(sam):
        if rec.is_paired:
            count += 1
        if i == num_reads:
            break
    if format_ratio(count, num_reads) > 0.75:
        return True
    elif format_ratio(count, num_reads) < 1 - paired_cutoff:
        return False
    else:
        raise RuntimeError("Unable to infer pairing from bamfile {}".format(path))


def index_fasta(fasta):
    """if fasta index exists, return path, otherwise generate and return path"""
    if not os.path.exists(fasta + '.fai'):
        cmd = ['samtools', 'faidx', fasta]
        runProc(cmd)
    return fasta + '.fai'


def validate_bam_fasta_pairs(bam_paths, fasta_index):
    """
    Make sure that this BAM is actually aligned to this fasta
    """
    fasta_sequences = {x.split()[0] for x in open(fasta_index)}
    for bam in bam_paths:
        handle = pysam.Samfile(bam)
        if set(handle.references) != fasta_sequences:
            base_err = 'Error: BAM {} does not have the same sequences as the genome FASTA {}'
            err = base_err.format(bam, fasta_index.replace('.fai', ''))
            raise RuntimeError(err)


def group_references(sam_handle, num_bases=10 ** 6, max_seqs=100):
    """
    Group up references by num_bases, unless that exceeds max_seqs.
    """
    name_iter = itertools.izip(*[sam_handle.references, sam_handle.lengths])
    name, size = name_iter.next()
    this_bin = [name]
    bin_base_count = size
    num_seqs = 1
    for name, size in name_iter:
        bin_base_count += size
        num_seqs += 1
        if bin_base_count >= num_bases or num_seqs > max_seqs:
            yield this_bin
            this_bin = [name]
            bin_base_count = size
            num_seqs = 1
        else:
            this_bin.append(name)


def main_hints_fn(target, bam_paths, db, genome, genome_fasta, hints_file):
    """
    Main driver function. Loops over each BAM, inferring paired-ness, then passing each BAM with one chromosome name
    for filtering. Maintains a chromosome based map of bams so that they can be merged prior to hint construction.
    """
    filtered_bam_tree = TempFileTree(get_tmp(target, global_dir=True, name="filter_file_tree"))
    fasta_index = index_fasta(genome_fasta)
    validate_bam_fasta_pairs(bam_paths, fasta_index)
    chrom_file_map = defaultdict(list)
    sam_handle = pysam.Samfile(bam_paths[0])
    references = list(group_references(sam_handle))
    for bam_path in bam_paths:
        paired = ["--paired", "--pairwiseAlignments"] if bam_is_paired(bam_path) is True else []
        for group in references:
            out_filter = filtered_bam_tree.getTempFile(suffix=".bam")
            chrom_file_map[tuple(group)].append(out_filter)
            target.addChildTargetFn(sort_by_name, args=[bam_path, group, out_filter, paired])
    target.setFollowOnTargetFn(merge_bam_wrapper, args=[chrom_file_map, genome, db, genome_fasta, hints_file])


def sort_by_name(target, bam_path, references, out_filter, paired):
    """
    Slices out a chromosome from a BAM, re-sorts by name, filters the reads, then re-sorts by position.
    filterBam does weird things when piped to stdout, so I don't do that.
    """
    tmp_filtered = get_tmp(target, name=".filtered.bam")
    cmd = [['samtools', 'view', '-b', bam_path],
           ['samtools', 'sort', '-O', 'bam', '-T', get_tmp(target), '-n', '-'],
           ['filterBam', '--uniq', '--in', '/dev/stdin', '--out', tmp_filtered]]
    cmd[-1].extend(paired)
    cmd[0].extend(references)
    runProc(cmd)
    cmd2 = ['samtools', 'sort', '-O', 'bam', '-T', get_tmp(target), tmp_filtered]
    runProc(cmd2, stdout=out_filter)
    runProc(['samtools', 'index', out_filter])


def merge_bam_wrapper(target, chrom_file_map, genome, db, genome_fasta, hints_file):
    """
    Wrapper to parallelize merging of bams by chromosome.
    """
    merged_bams = []
    for chrom, bam_files in chrom_file_map.iteritems():
        merged_path = get_tmp(target, global_dir=True, name='.combined.bam')
        target.addChildTargetFn(merge_bams, args=(bam_files, merged_path))
        merged_bams.append(merged_path)
    target.setFollowOnTargetFn(build_hints, args=[merged_bams, genome, db, genome_fasta, hints_file])


def merge_bams(target, bam_files, merged_path):
    fofn = get_tmp(target)
    with open(fofn, 'w') as outf:
        for x in bam_files:
            outf.write(x + "\n")
    cmd = ['samtools', 'merge', '-b', fofn, merged_path]
    runProc(cmd)


def build_hints(target, merged_bams, genome, db, genome_fasta, hints_file):
    """
    Driver function for hint building. Builts intron and exon hints, then calls cat_hints to do final concatenation
    and sorting.
    """
    if hints_file is None:
        hints_file = target.getGlobalTempDir()
    intron_hints_tree = TempFileTree(get_tmp(target, global_dir=True, name="intron_hints_tree"))
    exon_hints_tree = TempFileTree(get_tmp(target, global_dir=True, name="exon_hints_tree"))
    for merged_bam in merged_bams:
        intron_hints_path = intron_hints_tree.getTempFile(suffix=".intron.gff")
        target.addChildTargetFn(build_intron_hints, memory=8 * 1024 ** 3, args=[merged_bam, intron_hints_path])
        exon_hints_path = exon_hints_tree.getTempFile(suffix=".exon.gff")
        target.addChildTargetFn(build_exon_hints, memory=8 * 1024 ** 3, args=[merged_bam, exon_hints_path])
    target.setFollowOnTargetFn(cat_hints, args=[intron_hints_tree, exon_hints_tree, genome, db, genome_fasta,
                                                hints_file])


def build_exon_hints(target, bam_file, exon_gff_path):
    """
    Exon hints are built from peaks of coverage, using bam2wig and wig2hints.pl in the Augustus repo.
    """
    cmd = [['bam2wig', bam_file],
           ['wig2hints.pl', '--width=10', '--margin=10', '--minthresh=2', '--minscore=4', '--prune=0.1', '--src=W',
            '--type=ep', '--UCSC=/dev/null', '--radius=4.5', '--pri=4', '--strand="."']]
    runProc(cmd, stdout=exon_gff_path)


def build_intron_hints(target, bam_file, intron_hints_path):
    """
    Intron hints look for splice junctions in BAM files
    """
    cmd = ['bam2hints', '--intronsonly', '--in', bam_file, '--out', intron_hints_path]
    runProc(cmd)


def cat_hints(target, intron_hints_tree, exon_hints_tree, genome, db, genome_fasta, hints_file):
    """
    All intron and exon hint gff files are concatenated and then sorted.
    """
    if os.path.dirname(hints_file) != '':
        ensureDir(os.path.dirname(hints_file))
    all_gffs = intron_hints_tree.listFiles() + exon_hints_tree.listFiles()
    concat_hints = get_tmp(target, name="concat_hints")
    with open(concat_hints, 'w') as outf:
        for f in all_gffs:
            for l in open(f):
                outf.write(l)
    hints_tmp = tmpFileGet()
    # TODO: this takes forever. Surely this can be merged into one better sort command
    cmd = [['sort', '-n', '-k4,4', concat_hints],
           ['sort', '-s', '-n', '-k5,5'],
           ['sort', '-s', '-n', '-k3,3'],
           ['sort', '-s', '-k1,1'],
           ['join_mult_hints.pl']]
    runProc(cmd, stdout=hints_tmp)
    atomicInstall(hints_tmp, hints_file)
    target.setFollowOnTargetFn(load_db, args=[hints_file, db, genome, genome_fasta])


def load_db(target, hints_file, db, genome, genome_fasta, timeout=6000, intervals=30):
    """
    Final database loading.
    NOTE: Once done on all genomes, you want to run load2sqlitedb --makeIdx --dbaccess ${db}
    TODO: this times out sometimes. Rewrite to not be recursive.
    """
    def handle_concurrency(cmd, timeout, intervals, start_time=None):
        if start_time is None:
            start_time = time.time()
        elif time.time() - start_time >= timeout:
            raise RuntimeError("hints database still locked after {} seconds".format(timeout))
        p = subprocess.Popen(cmd, bufsize=-1, stderr=subprocess.PIPE)
        _, ret = p.communicate()
        if p.returncode == 0:
            return 1
        elif p.returncode == 1 and "locked" in ret:
            time.sleep(intervals)
            handle_concurrency(cmd, timeout, intervals, start_time)
        else:
            raise RuntimeError(ret)
    ensureFileDir(db)
    for f in [genome_fasta, hints_file]:
        cmd = ['load2sqlitedb', '--noIdx', '--species={}'.format(genome), '--dbaccess={}'.format(db), f]
        ret = handle_concurrency(cmd, timeout, intervals)


def external_main(args):
    assert os.path.exists(args.fasta)
    assert all([os.path.exists(x) for x in args.bams])
    args.defaultMemory = 8 * 1024 ** 3
    s = Stack(Target.makeTargetFn(main_hints_fn, args=[args.bams, args.database, args.genome, args.fasta,
                                                       args.hintsFile]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--hintsFile", default=None)
    parser.add_argument("--fasta", required=True)
    bamfiles = parser.add_mutually_exclusive_group(required=True)
    bamfiles.add_argument("--bamFiles", nargs="+", help="bamfiles being used", dest="bams")
    bamfiles.add_argument("--bamFofn", help="File containing list of bamfiles", dest="bams")
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    # below is an ugly hack to filter out tissues/centers by checking for the words in the file path
    assert os.path.exists(args.fasta)
    if not isinstance(args.bams, list):
        if not os.path.exists(args.bams):
            raise RuntimeError("ERROR: bamFofn does not exist.")
        bams = {x.rstrip() for x in open(args.bams)}
        args.bams = bams
    else:
        args.bams = set(args.bams)
    for bam in args.bams:
        assert os.path.exists(bam + 'bai'), 'Bam {} not indexed'.format(bam)
    args.defaultMemory = 8 * 1024 ** 3
    s = Stack(Target.makeTargetFn(main_hints_fn, args=[args.bams, args.database, args.genome, args.fasta,
                                                       args.hintsFile]))
    i = s.startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == '__main__':
    from comparativeAnnotator.augustus.build_hints_db import *
    main()
