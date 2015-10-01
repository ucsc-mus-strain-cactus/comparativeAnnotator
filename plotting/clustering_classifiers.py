# analyzing the classifiers failed most commonly in the basic set, then comparing this to the comprehensive set

import os
import sys
import itertools
import sqlite3 as sql
from lib.psl_lib import remove_alignment_number, remove_augustus_alignment_number
from lib.sqlite_lib import attach_database
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import seaborn as sns

def attach_databases(comp_ann_path):
    """
    Attaches all of the databases. Expects comp_ann_path to be the path that comparativeAnnotator wrote to.
    """
    con = sql.connect(os.path.join(comp_ann_path, "classify.db"))
    cur = con.cursor()
    attach_database(con, os.path.join(comp_ann_path, "augustusClassify.db"), "augustus")
    attach_database(con, os.path.join(comp_ann_path, "attributes.db"), "attributes")
    return con, cur

comp_ann_dir = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1504/comparativeAnnotation/2015-08-10_Augustus"
con, cur = attach_databases(comp_ann_dir)

tm_fields = ["CodingInsertions", "CodingDeletions", "StartOutOfFrame", "FrameShift",
                      "AlignmentAbutsLeft", "AlignmentAbutsRight", "AlignmentPartialMap", "BadFrame", "BeginStart",
                      "CdsGap", "CdsMult3Gap", "UtrGap", "UnknownGap", "CdsUnknownSplice", "UtrUnknownSplice",
                      "EndStop", "InFrameStop", "ShortCds", "UnknownBases", "AlignmentAbutsUnknownBases"]


aug_fields = ['AugustusParalogy', 'AugustusExonGain', 'AugustusExonLoss', 'AugustusNotSameStrand',
              'AugustusNotSameStartStop', 'AugustusNotSimilarTerminalExonBoundaries',
              'AugustusNotSimilarInternalExonBoundaries']


tm_fail_counts = {}
for tm_classifier in tm_fields:
    cmd = """SELECT main.'C57B6NJ'.'AlignmentId' FROM main.'C57B6NJ' WHERE main.'C57B6NJ'.'{}' = 1""".format(tm_classifier)
    r = cur.execute(cmd).fetchall()
    tm_fail_counts[tm_classifier] = len(r)


aug_fail_counts = {}
for aug_classifier in aug_fields:
    cmd = """SELECT augustus.'C57B6NJ'.'AlignmentId' FROM augustus.'C57B6NJ' WHERE augustus.'C57B6NJ'.'{}' = 1""".format(aug_classifier)
    r = cur.execute(cmd).fetchall()
    aug_fail_counts[aug_classifier] = len(r)




# decided to do heiarchaical clustering in R instead of dealing with matplotlib

#sqlite3 command to dump to text: select AlignmentId,CodingInsertions,CodingDeletions,StartOutOfFrame,FrameShift,AlignmentAbutsLeft,AlignmentAbutsRight,AlignmentPartialMap,BadFrame,BeginStart,CdsGap,CdsMult3Gap,UtrGap,UnknownGap,CdsUnknownSplice,UtrUnknownSplice,EndStop,InFrameStop,ShortCds,UnknownBases,AlignmentAbutsUnknownBases from C57B6NJ;

# R code

data <- read.csv("~/mus_strain_data/pipeline_data/comparative/1504/comparativeAnnotation/2015-08-10_Augustus/C57B6NJ.transMap.csv", row.names=1, header=T, na.strings="")
# need to convert NA to 0
data[is.na(data)] <- 0
mat <- sapply(as.data.frame(data), as.logical)
mat.t <- t(mat)
library(stats)
d <- dist(mat.t, method="binary")
hc <- hclust(d, method="ward")
pdf("transMap_clustered_classifiers.pdf")
plot(hc)
dev.off()

pdf("transMap_categories.pdf")
data.transactions <- as(mat, "transactions")
itemFrequencyPlot(data.transactions, support=0.1, cex.names=0.7)
dev.off()

rules <- apriori(data.transactions, parameter=list(support=0.01, confidence=0.7))
# i have no fucking clue how to interpret this, but at least the frequency plot is nice


# now I want to re-run this on just Basic set protein_coding. I can't just dump to csv this time.
from scripts.coverage_identity_ok_plots import *
from scripts.consensus import *
attrs = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeAttrsVM4.tsv"
coding_ids = get_all_ids(attrs, biotype="protein_coding")
basic_ids = {x.split()[0] for x in open("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeBasicVM4.gp")}
basic_coding = {x for x in basic_ids if x in coding_ids}
tm_cmd = """SELECT AlignmentId,{} FROM main.'C57B6NJ'""".format(",".join(tm_fields))
r = cur.execute(tm_cmd).fetchall()
r_coding = [x for x in r if remove_alignment_number(x[0]) in basic_coding]
with open("transmap_coding_only.csv", "w") as outf:
    outf.write("AlignmentId," + ",".join(tm_fields) + "\n")
    for x in r_coding:
        outf.write(",".join(map(str, x)) + "\n")


data <- read.csv("/hive/users/ifiddes/comparativeAnnotator/transmap_coding_only.csv", row.names=1, header=T, na.strings="None")
data[is.na(data)] <- 0
mat <- sapply(as.data.frame(data), as.logical)
mat.t <- t(mat)
library(stats)
d <- dist(mat.t, method="binary")
hc <- hclust(d, method="ward")
pdf("transMap_clustered_classifiers_coding_basic.pdf")
plot(hc)
dev.off()
pdf("transMap_categories_coding_basic.pdf")
data.transactions <- as(mat, "transactions")
itemFrequencyPlot(data.transactions, support=0.1, cex.names=0.7)
dev.off()


data <- read.csv("~/mus_strain_data/pipeline_data/comparative/1504/comparativeAnnotation/2015-08-10_Augustus/C57B6NJ.augustus.csv", row.names=1, header=T, na.strings="")
# need to convert NA to 0
data[is.na(data)] <- 0
mat <- sapply(as.data.frame(data), as.logical)
mat.t <- t(mat)
library(stats)
d <- dist(mat.t, method="binary")
hc <- hclust(d, method="ward")
pdf("augustus_clustered_classifiers.pdf")
plot(hc)
dev.off()

pdf("augustus_categories.pdf")
data.transactions <- as(mat, "transactions")
itemFrequencyPlot(data.transactions, support=0.1, cex.names=0.7)
dev.off()


# now I want the complement - transcripts in comprehensive that are not in basic
from scripts.coverage_identity_ok_plots import *
from scripts.consensus import *
attrs = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeAttrsVM4.tsv"
coding_ids = get_all_ids(attrs, biotype="protein_coding")
basic_ids = {x.split()[0] for x in open("/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeBasicVM4.gp")}
basic_coding = {x for x in basic_ids if x in coding_ids}
comp_specific_ids = coding_ids - basic_ids
tm_cmd = """SELECT AlignmentId,{} FROM main.'C57B6NJ'""".format(",".join(tm_fields))
r = cur.execute(tm_cmd).fetchall()
r_coding = [x for x in r if remove_alignment_number(x[0]) in comp_specific_ids]
with open("transmap_comp_specific.csv", "w") as outf:
    outf.write("AlignmentId," + ",".join(tm_fields) + "\n")
    for x in r_coding:
        outf.write(",".join(map(str, x)) + "\n")


data <- read.csv("/hive/users/ifiddes/comparativeAnnotator/transmap_comp_specific.csv", row.names=1, header=T, na.strings="None")
# need to convert NA to 0
data[is.na(data)] <- 0
mat <- sapply(as.data.frame(data), as.logical)
mat.t <- t(mat)
library(stats)
d <- dist(mat.t, method="binary")
hc <- hclust(d, method="ward")
pdf("transMap_comp_only_classifiers.pdf")
plot(hc)
dev.off()

pdf("transMap_comp_only_categories.pdf")
data.transactions <- as(mat, "transactions")
itemFrequencyPlot(data.transactions, support=0.1, cex.names=0.7)
dev.off()

# not sure what I want to use this for
tm_grouped_fields = [["CodingInsertions", "CodingDeletions"], ["CodingMult3Deletions", "CodingMult3Insertions"], "StartOutOfFrame", "FrameShift",
                     ["AlignmentAbutsLeft", "AlignmentAbutsRight"], "AlignmentPartialMap", "BadFrame", ["BeginStart", "EndStop"],
                     ["CdsGap", "CdsMult3Gap"], "UtrGap", "UnknownGap", ["CdsUnknownSplice", "UtrUnknownSplice"], "InFrameStop", "ShortCds",
                     ["UnknownBases", "AlignmentAbutsUnknownBases"]]


# how many genes are present in comp that are not present in basic?
attr_map = {x.split()[3]: x.split()[0] for x in open(attrs)}
comp_genes = {attr_map[x] for x in attr_map if x not in basic_ids}
# 19,909. What biotypes do we have here?
biotype_map = {x.split()[0]: x.split()[2] for x in open(attrs)}
from collections import Counter
biotypes = Counter()
for g in comp_genes:
    biotypes[biotype_map[g]] += 1
# Counter({'protein_coding': 10545, 'processed_pseudogene': 5421, 'unprocessed_pseudogene': 2042, 'lincRNA': 829, 'processed_transcript': 370, 'pseudogene': 179, 'transcribed_processed_pseudogene': 139, 'transcribed_unprocessed_pseudogene': 128, 'IG_V_pseudogene': 69, 'IG_V_gene': 40, 'TEC': 33, 'TR_V_gene': 23, 'antisense': 21, 'unitary_pseudogene': 15, 'TR_V_pseudogene': 14, 'polymorphic_pseudogene': 9, 'sense_intronic': 8, 'IG_C_gene': 6, 'IG_D_pseudogene': 4, 'sense_overlapping': 3, 'TR_J_gene': 2, 'IG_C_pseudogene': 1, 'translated_processed_pseudogene': 1, 'geneType': 1, 'IG_J_gene': 1, 'TR_J_pseudogene': 1, '3prime_overlapping_ncrna': 1, 'translated_unprocessed_pseudogene': 1, 'IG_D_gene': 1, 'TR_C_gene': 1})
# still over half coding, but highly enriched for pseudogenes as expected

