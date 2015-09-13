import sys
sys.path.append("/hive/users/ifiddes/comparativeAnnotator")
import os
from collections import defaultdict
from lib.sequence_lib import GenePredTranscript, getGenePredTranscripts, transcriptListToDict
from scripts.plot_functions import skip_header

attr_path = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeAttrsVM4.tsv"
gp_path = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/comparativeAnnotation/2015-09-01_Augustus/consensus/geneSets/129S1consensusGeneSet.gp"


def parse_attributes(attr_path):
    tx_gene_map = defaultdict(list)
    ens_gene_map = {}
    for x in skip_header(attr_path):
        gene_id, gene_name, gene_type, transcript_id, transcript_type = x.split()
        tx_gene_map[gene_id].append(transcript_id)
        ens_gene_map[gene_id] = gene_name
    return tx_gene_map, ens_gene_map



transcripts = getGenePredTranscripts(gp_path)
transcript_dict = transcriptListToDict(transcripts, noDuplicates=True)
tx_gene_map, ens_gene_map = parse_attributes(attr_path)

