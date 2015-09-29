from scripts.plot_functions import *
from scripts.consensus import *


comp_ann_path = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/comparativeAnnotation/2015-09-01_Augustus"
attr_path = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeAttrsVM4.tsv"
genomes = "129S1 AJ AKRJ BALBcJ C3HHeJ C57B6NJ CASTEiJ CBAJ DBA2J FVBNJ LPJ NODShiLtJ NZOHlLtJ PWKPhJ SPRETEiJ WSBEiJ CAROLIEiJ PAHARIEiJ".split()
tm_gp_base = "/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/comparative/1504/transMap/2015-05-28/transMap/{}/syn/transMapGencodeCompVM4.gp"
aug_gp_base = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/augustus/tmr/{}.gp"
basic_gp = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeBasicVM4.gp"
comp_gp = "/cluster/home/ifiddes/mus_strain_data/pipeline_data/comparative/1504/transMap/2015-05-28/data/wgEncodeGencodeCompVM4.gp"
out_path = "."
biotype = "protein_coding"


gene_map = get_gene_map(attr_path)
gene_biotype_map = get_gene_biotype_map(attr_path)


stats_dir = os.path.join(comp_ann_path, "augustus_stats")
con, cur = attach_databases(comp_ann_path, has_augustus=True)

binned_transcript_holder = {}
for genome in genomes:
    coding_ok = get_all_ok(cur, genome, tm_coding_classifiers)
    tm_gp = tm_gp_base.format(genome)
    aug_gp = aug_gp_base.format(genome)
    gps = load_gps([tm_gp, aug_gp])
    chr_y_ids = gp_chrom_filter(comp_gp)
    ens_ids = get_all_ids(attr_path, biotype=biotype) - chr_y_ids
    reverse_name_map = get_reverse_name_map(cur, genome, whitelist=ens_ids, has_augustus=True)
    stats_dict = merge_stats(cur, stats_dir, genome)
    binned_transcripts = bin_transcripts(reverse_name_map, stats_dict, coding_ok, ens_ids)
    binned_transcript_holder[genome] = binned_transcripts
    break


comp_ids = get_gp_ids(comp_gp)
basic_ids = get_gp_ids(basic_gp)
coding_ids = get_all_ids(attr_path, biotype="protein_coding")
basic_coding = basic_ids & coding_ids
comp_coding = comp_ids & coding_ids
complement_coding = comp_coding - basic_coding
genome_order = hard_coded_genome_order
cat = "Basic"
ids = basic_coding
title_string = base_title_string.format(len(ids), title_string_dict[cat].format("Transcript"))
out_name = file_name_dict[cat]
coding_metrics = OrderedDict()
for g in genome_order:
    bins = binned_transcript_holder[g]
    coding_metrics[g] = make_counts_frequency(make_tx_counts_dict(bins, filter_set=ens_ids))


categories = zip(*coding_metrics[g])[0]
results = [[g, zip(*coding_metrics[g])[1]] for g in coding_metrics]

stacked_barplot(results, categories, out_path, out_name, title_string, color_palette=paired_palette)

binned_transcripts = {"bestOk": [], "augAltOk": [], "tmAltOk": [], "bestNotOk": [], "augAltNotOk": [], 
                      "tmAltNotOk": [], "fail": [], "discarded": [], "tieIds": set()}


for ens_id in ens_ids:
    aln_ids = reverse_name_map[ens_id]
    ok_candidates, not_ok_candidates, discarded = find_ok_not_ok_candidates(stats_dict, aln_ids, ok_ids,
                                                                            discard_cov_cutoff=0.50,
                                                                            filter_cov_cutoff=0.80)
    if len(ok_candidates) > 0:
        ok_id, aug_alt, tm_alt, is_tie = analyze_candidates(ok_candidates)
        if is_tie:
            break



def analyze_candidates(candidates, aug_re=re.compile("I[1-2]-")):
    """
    Analyzes candidate transcripts, finding the winner and splitting the alternatives based on coming from augustus
    or transMap. If there are multiple equal winners from both transMap and Augustus, reports so.
    aug_re is used if Augustus was ran in two modes (I1 and I2) to prevent counting that as a tie with transMap.
    """
    winner_ids = find_best_aln(candidates)
    winner_id = winner_ids[0]
    if len({"-".join(aug_re.split(x[0])) for x in winner_ids}) > 1:
        is_tie = True
    else:
        is_tie = False
    aug_alts, tm_alts = split_alternatives(candidates, winner_id)
    return winner_id, aug_alts, tm_alts, is_tie
