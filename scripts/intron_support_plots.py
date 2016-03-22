genomes = ['NZO_HlLtJ', 'C57BL_6NJ', 'LP_J', 'C3H_HeJ', 'FVB_NJ', 'CBA_J', 'CAROLI_EiJ', 'PWK_PhJ', 'SPRET_EiJ',  'A_J', '129S1_SvImJ', 'Pahari_EiJ', 'NOD_ShiLtJ', 'AKR_J', 'WSB_EiJ', 'DBA_2J', 'CAST_EiJ', 'BALB_cJ']

for genome in genomes:
    bed_dir = '/hive/groups/recon/projs/mus_strain_cactus/pipeline_data/rnaseq/munged_STAR_data/REL-1509-chromosomes/{}/'.format(genome)
    beds = [os.path.join(bed_dir, x) for x in listdir(bed_dir, suffix='.sj.bed')]
    assert len(beds) > 0
    cmd = ['python', 'submodules/comparativeAnnotator/scripts/validate_splice_junctions.py', '--consensus_gp',
    'mouse_output/CGP_consensus/{}.CGP_consensus.gp'.format(genome), '--star_junctions']
    cmd.extend(beds)
    cmd.extend(['--out', 'supported_data/{}.txt'.format(genome)])
    runProc(cmd)


import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


for genome in genomes:
    df = pd.read_csv('{}.txt'.format(genome), header=0, sep='\t')
    p = sns.pairplot(df, x_vars=['NumIntrons'], y_vars=['NumSupported'], hue='Source', size=10, kind='reg')
    ax = p.axes[0][0]
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 120)
    plt.title("Splice support in {}".format(genome))
    plt.savefig('{}_pair.png'.format(genome), format='png')
    plt.close('all')
    p2 = sns.distplot(df.Ratio.dropna(), norm_hist=False, kde=False, bins=20)
    p2.set_xlabel('Percentage of splices supported by RNAseq')
    p2.set_ylabel('Number of transcripts')
    plt.title("Ratio of splice support in {}".format(genome))
    plt.savefig('{}_hist.png'.format(genome), format='png')
    plt.close('all')
    p3 = sns.pairplot(df, x_vars=['NumIntrons'], y_vars=['NumSupported'], size=10, plot_kws={'alpha':0.3})
    p3.map(sns.kdeplot, cmap="Blues_d", n_levels=50)
    ax = p3.axes[0][0]
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 120)
    plt.title("Splice support in {}".format(genome))
    plt.savefig('{}_pair_kde.png'.format(genome), format='png')
    plt.close('all')
    df_20 = df[df.NumIntrons <= 20]
    p4 = sns.pairplot(df_20, x_vars=['NumIntrons'], y_vars=['NumSupported'], size=10, plot_kws={'alpha':0.3})
    p4.map(sns.kdeplot, cmap="Blues_d", n_levels=50)
    ax = p4.axes[0][0]
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 20)
    plt.title("Splice support in {}".format(genome))
    plt.savefig('{}_pair_kde_zoomed.png'.format(genome), format='png')
    plt.close('all')
    df_50 = df[df.NumIntrons <= 50]
    p5 = sns.distplot(df_50.NumIntrons.dropna(), norm_hist=False, kde=False, bins=20)
    plt.title('Number of introns per transcript in {}'.format(genome))
    p5.set_xlabel('Number of introns')
    p5.set_ylabel('Number of transcripts')
    plt.savefig('{}_num_introns.png'.format(genome), format='png')
    plt.close('all')


genome = 'C57B6J'
ref_df = pd.read_csv('{}.txt'.format(genome), header=0, sep='\t')
p = sns.pairplot(df, x_vars=['NumIntrons'], y_vars=['NumSupported'], size=10, kind='reg')
ax = p.axes[0][0]
ax.set_xlim(0, 120)
ax.set_ylim(0, 120)
plt.title("Splice support in {}".format(genome))
plt.savefig('{}_pair.png'.format(genome), format='png')
plt.close('all')
p2 = sns.distplot(ref_df.Ratio.dropna(), norm_hist=False, kde=False, bins=20)
p2.set_xlabel('Percentage of splices supported by RNAseq')
p2.set_ylabel("Number of transcripts")
plt.title("Ratio of splice support in {}".format(genome))
plt.savefig('{}_hist.png'.format(genome), format='png')
plt.close('all')
p3 = sns.pairplot(ref_df, x_vars=['NumIntrons'], y_vars=['NumSupported'], size=10, plot_kws={'alpha':0.3})
p3.map(sns.kdeplot, cmap="Blues_d", n_levels=50)
ax = p3.axes[0][0]
ax.set_xlim(0, 120)
ax.set_ylim(0, 120)
plt.title("Splice support in {}".format(genome))
plt.savefig('{}_pair_kde.png'.format(genome), format='png')
plt.close('all')
ref_df_20 = ref_df[ref_df.NumIntrons <= 20]
p4 = sns.pairplot(ref_df_20, x_vars=['NumIntrons'], y_vars=['NumSupported'], size=10, plot_kws={'alpha':0.3})
p4.map(sns.kdeplot, cmap="Blues_d", n_levels=50)
ax = p4.axes[0][0]
ax.set_xlim(0, 20)
ax.set_ylim(0, 20)
plt.title("Splice support in {}".format(genome))
plt.savefig('{}_pair_kde_zoomed.png'.format(genome), format='png')
plt.close('all')
ref_df_50 = ref_df[ref_df.NumIntrons <= 50]
p5 = sns.distplot(ref_df_50.NumIntrons.dropna(), norm_hist=False, kde=False, bins=20)
plt.title('Number of introns per transcript in {}'.format(genome))
p5.set_xlabel('Number of introns')
p5.set_ylabel('Number of transcripts')
plt.savefig('{}_num_introns.png'.format(genome), format='png')
plt.close('all')


for genome in genomes:
    df = pd.read_csv('{}.txt'.format(genome), header=0, sep='\t')
    c_df = pd.merge(df, ref_df, how='left', on='TxId', suffixes=("", "_ref"))
    p = sns.pairplot(c_df, x_vars=['NumIntrons_ref'], y_vars=['NumIntrons'], size=10, kind='reg')
    ax = p.axes[0][0]
    ax.set_xlim(0, 120)
    ax.set_ylim(0, 120)
    plt.title("Change in number of introns from reference in {}".format(genome))
    plt.savefig('{}_change_intron_number.png'.format(genome), format='png')
    plt.close('all')


for genome in genomes:
    files = [x for x in os.listdir('.') if genome in x]
    os.mkdir(genome)
    for f in files:
        os.rename(f, os.path.join(genome, f))
