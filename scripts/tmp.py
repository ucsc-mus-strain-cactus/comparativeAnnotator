from scripts.plot_functions import *
from scripts.consensus import *
import cPickle as pickle
binned_transcript_holder = pickle.load(open("bins.pickle"))
comp_ann_path = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/comparativeAnnotation/2015-09-01_Augustus"
attr_path = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeAttrsVM4.tsv"
genomes = "129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ CAROLIEiJ PAHARIEiJ".split()
tm_gp_base = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1504/transMap/2015-05-28/transMap/{}/transMapGencodeCompVM4.gp"
aug_gp_base = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/augustus/tmr/{}.gp"
basic_gp = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeBasicVM4.gp"
comp_gp = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeCompVM4.gp"
out_path = "."
biotype = "protein_coding"
genome_order = make_coding_transcript_plots(binned_transcript_holder, out_path, comp_gp, basic_gp, attr_path)
gene_map = get_gene_map(attr_path)
ok_gene_by_biotype(binned_transcript_holder, out_path, attr_path, gene_map, genome_order, biotype)



out_path = os.path.join(comp_ann_path, "consensus")
stats_dir = os.path.join(comp_ann_path, "augustus_stats")
con, cur = attach_databases(comp_ann_path)
biotypes = get_all_biotypes(attr_path)
gene_map = get_gene_map(attr_path)
gene_biotype_map = get_gene_biotype_map(attr_path)

binned_transcript_holder = defaultdict(dict)
for genome in genomes:
    coding_ok = get_all_ok(cur, genome, tm_coding_classifiers)
    noncoding_ok = get_all_ok(cur, genome, tm_noncoding_classifiers)
    tm_gp = tm_gp_base.format(genome)
    aug_gp = aug_gp_base.format(genome)
    gps = load_gps([tm_gp, aug_gp])
    #biotypes = ["protein_coding"]
    biotypes = ["lincRNA", "miRNA", "snoRNA"]
    for biotype in biotypes:
        ens_ids = get_all_ids(attr_path, biotype=biotype)
        reverse_name_map = get_reverse_name_map(cur, genome, ens_ids)
        stats_dict = merge_stats(cur, stats_dir, genome)
        binned_transcripts = bin_transcripts(reverse_name_map, stats_dict, coding_ok, ens_ids)
        binned_transcript_holder[genome][biotype] = binned_transcripts


