from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from lib.psl_lib import PslRow, removeAugustusAlignmentNumber, removeAlignmentNumber
from sonLib.bioio import fastaRead, fastaWrite, popenCatch, system, getRandomAlphaNumericString
from pyfaidx import Fasta
from lib.general_lib import formatRatio
import errno
import os
import argparse


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def coverage(p_list):
    m = sum(x.matches for x in p_list)
    mi = sum(x.misMatches for x in p_list)
    rep = sum(x.repMatches for x in p_list)
    return formatRatio(m + mi + rep, p_list[0].qSize)


def identity(p_list):
    m = sum(x.matches for x in p_list)
    mi = sum(x.misMatches for x in p_list)
    rep = sum(x.repMatches for x in p_list)
    ins = sum(x.qNumInsert for x in p_list)
    return formatRatio(m + rep, m + rep + mi + ins)


def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))


def align(target, g, target_fasta, chunk, ref_fasta, out_path):
    g_f = Fasta(target_fasta)
    r_f = Fasta(ref_fasta)
    results = []
    for aug_aId in chunk:
        aId = removeAugustusAlignmentNumber(aug_aId)
        gencode_id = removeAlignmentNumber(aId)
        gencode_seq = str(r_f[gencode_id])
        aug_seq = str(g_f[aug_aId])
        tmp_aug = os.path.join(target.getLocalTempDir(), "tmp_aug")
        tmp_gencode = os.path.join(target.getLocalTempDir(), "tmp_gencode")
        fastaWrite(tmp_aug, aug_aId, aug_seq)
        fastaWrite(tmp_gencode, gencode_id, gencode_seq)
        r = popenCatch("blat {} {} -out=psl -noHead /dev/stdout".format(tmp_gencode, tmp_aug))
        r = r.split("\n")[:-3]
        if len(r) == 0:
            results.append([aug_aId, "0", "0"])
        else:
            p_list = [PslRow(x) for x in r]
            results.append(map(str, [aug_aId, identity(p_list), coverage(p_list)]))
    with open(os.path.join(out_path, getRandomAlphaNumericString(10) + ".txt"), "w") as outf:
        for x in results:
            outf.write("\t".join(x) + "\n")


def cat(target, g, in_path, out_dir):
    in_p = os.path.join(in_path, "*")
    out_p = os.path.join(out_dir, "augustus", g + ".stats")
    system("cat {} > {}".format(in_p, out_p))


def wrapper(target, genomes, ref_fasta, out_dir):
    for g in genomes:
        out_path = os.path.join(out_dir, "tmp", g)
        mkdir_p(out_path)
        for f in os.listdir(out_path):
            os.remove(os.path.join(out_path, f))
        target_fasta = os.path.join(out_dir, "augustus", g + ".fa")
        faidx = os.path.join(out_dir, "augustus", g + ".fa.fai")
        aug_aIds = [x.split()[0] for x in open(faidx)]
        for chunk in chunker(aug_aIds, 200):
            target.addChildTargetFn(align, args=(g, target_fasta, chunk, ref_fasta, out_path))
    target.setFollowOnTargetFn(wrapper2, args=(genomes, out_dir))


def wrapper2(target, genomes, out_dir):
    for g in genomes:
        out_path = os.path.join(out_dir, "tmp", g)
        target.addChildTargetFn(cat, args=(g, out_path, out_dir))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--refFasta", required=True)
    parser.add_argument("--outDir", required=True)
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    i = Stack(Target.makeTargetFn(wrapper, args=(args.genomes, args.refFasta, 
                                                 args.outDir))).startJobTree(args)
    if i != 0:
        raise RuntimeError("Got failed jobs")
    system("rm -rf {}".format(os.path.join(args.outDir, "tmp")))


if __name__ == '__main__':
    from scripts.alignAugustus import *
    main()
