"""
Wrapper to run CGP consensus. In the future, will be integrated in the pipeline once CGP prediction construction
is also integrated.
"""
import sys
import os
import argparse
import luigi
import csv
from frozendict import frozendict
os.environ['PYTHONPATH'] = './:./submodules:./submodules/pycbio:./submodules/comparativeAnnotator'
sys.path.extend(['./', './submodules', './submodules/pycbio', './submodules/comparativeAnnotator'])
from jobTree.scriptTree.stack import Stack
from lib.parsing import HashableNamespace, FileArgumentParser
from pipeline.abstract_classes import AbstractJobTreeTask, AbstractAtomicFileTask, RowExistsSqlTarget
from comparativeAnnotator.scripts.align_cgp_cds import align_cgp_cds
from comparativeAnnotator.scripts.cgp_consensus import cgp_consensus
from comparativeAnnotator.scripts.cgp_consensus_plots import generate_consensus_plots


class CgpConsensusNamespace(HashableNamespace):
    """
    Add a repr to prevent spamming the luigi logfile
    """
    def __repr__(self):
        return 'CgpConsensus-{}'.format(self.genome)


class AlignCgpNamespace(HashableNamespace):
    """
    Add a repr to prevent spamming the luigi logfile
    """
    def __repr__(self):
        return 'AlignCgp-{}'.format(self.genome)


class AlignCdsNamespace(HashableNamespace):
    """
    Add a repr to prevent spamming the luigi logfile
    """
    def __repr__(self):
        return 'AlignCds-{}'.format(self.genome)


class CgpPlotNamespace(HashableNamespace):
    """
    Add a repr to prevent spamming the luigi logfile
    """
    def __repr__(self):
        return 'CgpConsensusPlots-{}'.format(','.join(self.genomes))


class CgpConsensus(luigi.WrapperTask):
    """
    Main driver class.
    """
    params = luigi.Parameter()

    def create_args(self):
        args_holder = {}
        for genome, (consensus_gp, cgp_gp, genome_fasta, cgp_intron_bits) in self.params.file_map.iteritems():
            args = CgpConsensusNamespace()
            args.genome = genome
            args.ref_genome = self.params.refGenome
            args.norestart = self.params.norestart
            args.comp_db = self.params.compDb
            args.cgp_db = os.path.join(self.params.workDir, 'comparativeAnnotator', 'cgp_stats.db')
            args.consensus_gp = consensus_gp
            args.cgp_gp = cgp_gp
            args.cgp_intron_bits = cgp_intron_bits
            args.metrics_dir = os.path.join(self.params.workDir, 'CGP_consensus_metrics')
            # Alignment shared args
            tmp = HashableNamespace(**vars(self.params.jobTreeOptions))
            tmp.refGenome = self.params.refGenome
            tmp.genome = genome
            tmp.refTranscriptFasta = self.params.refTranscriptFasta
            tmp.compDb = self.params.compDb
            tmp.cgpDb = args.cgp_db
            tmp.targetGenomeFasta = genome_fasta
            tmp.defaultMemory = 8 * 1024 ** 3
            # align CGP
            args.align_cgp = AlignCgpNamespace(**vars(tmp))
            args.align_cgp.jobTree = os.path.join(self.params.jobTreeDir, 'alignCgp', genome)
            args.align_cgp.mode = 'cgp'
            args.align_cgp.gp = cgp_gp
            args.align_cgp.table = '_'.join([genome, args.align_cgp.mode])
            # align consensus CDS
            args.align_cds = AlignCdsNamespace(**vars(tmp))
            args.align_cds.jobTree = os.path.join(self.params.jobTreeDir, 'alignCds', genome)
            args.align_cds.mode = 'consensus'
            args.align_cds.gp = consensus_gp
            args.align_cds.table = '_'.join([genome, args.align_cds.mode])
            args.align_cds.genome = args.align_cgp.genome = genome
            # output
            args.output_gp = os.path.join(self.params.outputDir, 'CGP_consensus', genome + '.CGP_consensus.gp')
            args.output_gtf = os.path.join(self.params.outputDir, 'CGP_consensus', genome + '.CGP_consensus.gtf')
            args_holder[genome] = args
        return args_holder

    def create_plot_args(self):
        plot_args = CgpPlotNamespace()
        args_holder = self.create_args()
        plot_args.args_holder = frozendict(args_holder)
        plot_args.genomes = self.params.targetGenomes
        plot_args.metrics_dir = os.path.join(self.params.workDir, 'CGP_consensus_metrics')
        plot_args.plot_dir = os.path.join(self.params.outputDir, 'CGP_consensus_plots')
        plot_args.addition_plot = os.path.join(plot_args.plot_dir, 'cgp_addition.pdf')
        plot_args.replace_plot = os.path.join(plot_args.plot_dir, 'cgp_replace_rate.pdf')
        plot_args.new_isoform_plot = os.path.join(plot_args.plot_dir, 'new_isoform_rate.pdf')
        plot_args.missing_plot = os.path.join(plot_args.plot_dir, 'missing_genes_rescued.pdf')
        plot_args.consensus_gene_plot = os.path.join(plot_args.plot_dir, 'consensus_gene_plot.pdf')
        plot_args.consensus_tx_plot = os.path.join(plot_args.plot_dir, 'consensus_tx_plot.pdf')
        plot_args.consensus_tx_plot = os.path.join(plot_args.plot_dir, 'consensus_cgp_match_removed.pdf')
        plot_args.plots = (plot_args.addition_plot, plot_args.replace_plot, plot_args.new_isoform_plot,
                           plot_args.missing_plot, plot_args.consensus_gene_plot, plot_args.consensus_tx_plot,
                           plot_args.consensus_tx_plot)
        return plot_args, args_holder

    def requires(self):
        plot_args, args_holder = self.create_plot_args()
        for args in args_holder.itervalues():
            yield Align(args.align_cgp)
            yield Align(args.align_cds)
            yield GenerateConsensus(args, target_file=args.output_gp)
            yield ConvertGpToGtf(args, target_file=args.output_gtf)
        yield CgpConsensusPlots(plot_args)


class Align(AbstractJobTreeTask):
    def output(self):
        row_query = 'SELECT aln_table FROM completionFlags WHERE aln_table = "{}"'.format(self.cfg.table)
        return RowExistsSqlTarget(self.cfg.cgpDb, 'completionFlags', row_query)

    def run(self):
        self.start_jobtree(self.cfg, align_cgp_cds)


class GenerateConsensus(AbstractAtomicFileTask):
    def requires(self):
        return Align(cfg=self.cfg.align_cgp), Align(cfg=self.cfg.align_cds)

    def run(self):
        consensus_records = cgp_consensus(self.cfg)
        out_h = self.output().open('w')
        for rec in consensus_records:
            out_h.write(rec)
        out_h.close()


class ConvertGpToGtf(AbstractAtomicFileTask):
    """
    Converts the output genePred to GTF
    """
    def requires(self):
        return GenerateConsensus(cfg=self.cfg, target_file=self.cfg.output_gp)

    def run(self):
        cmd = [['bin/fixGenePredScore', self.cfg.output_gp],
               ['genePredToGtf', '-honorCdsStat', '-utr', 'file', '/dev/stdin', '/dev/stdout']]
        self.run_cmd(cmd)


class CgpConsensusPlots(luigi.Task):
    plot_args = luigi.Parameter()

    def requires(self):
        r = []
        for cfg in self.plot_args.args_holder.values():
            r.append(GenerateConsensus(cfg=cfg, target_file=cfg.align_cgp.gp))
            r.append(GenerateConsensus(cfg=cfg, target_file=cfg.align_cds.gp))
        return r

    def output(self):
        return [luigi.LocalTarget(x) for x in self.plot_args.plots]

    def run(self):
        generate_consensus_plots(self.plot_args)


def parse_args():
    """
    Build argparse object, parse arguments. See the parsing library for a lot of the features used here.
    """
    parser = FileArgumentParser(description=__doc__)
    parser.add_argument_with_check('--config', metavar='FILE', help='TSV file containing all inputs. Each line should '
                                                                    'have the columns genome, consensusGenePred, '
                                                                    'cgpGenePred, genomeFasta. cgpIntronBits is an '
                                                                    'optional field if you have it.')
    parser.add_argument('--refGenome', required=True, help='Reference genome')
    parser.add_argument('--refTranscriptFasta', required=True, help='Reference gene set transcript FASTA')
    parser.add_argument('--targetGenomes', nargs='+', required=True, help='Ordered target genomes.')
    parser.add_argument_with_check('--compDb', required=True, metavar='FILE', help='comparativeAnnotator database.')
    parser.add_argument_with_mkdir_p('--jobTreeDir', default='jobTrees', metavar='DIR',
                                     help='Work directory. Will contain intermediate files that may be useful.')
    parser.add_argument_with_mkdir_p('--outputDir', default='output', metavar='DIR',
                                     help='Output directory.')
    parser.add_argument_with_mkdir_p('--workDir', default='work', metavar='DIR',
                                     help='Work directory.')
    parser.add_argument('--norestart', action='store_true', default=False,
                        help='Set to force jobtree pipeline components to start from the beginning instead of '
                             'attempting a restart.')
    parser.add_argument('--localCores', default=12, metavar='INT',
                        help='Number of local cores to use. (default: %(default)d)')
    jobtree = parser.add_argument_group('jobTree options. Read the jobTree documentation for other options not shown')
    jobtree.add_argument('--batchSystem', default='parasol', help='jobTree batch system.')
    jobtree.add_argument('--parasolCommand', default='./bin/remparasol',
                         help='Parasol command used by jobTree. Used to remap to host node.')
    jobtree.add_argument('--maxThreads', default=4,
                         help='maxThreads for jobTree. If not using a cluster, this should be --localCores/# genomes')
    jobtree_parser = argparse.ArgumentParser(add_help=False)
    Stack.addJobTreeOptions(jobtree_parser)
    args = parser.parse_args(namespace=HashableNamespace())
    args.targetGenomes = tuple(args.targetGenomes)
    args.jobTreeOptions = jobtree_parser.parse_known_args(namespace=HashableNamespace())[0]
    args.jobTreeOptions.jobTree = None
    args.jobTreeOptions.__dict__.update({x: y for x, y in vars(args).iteritems() if x in args.jobTreeOptions})
    # munge parsed args, verify, make hashable
    file_map = {}
    for rec in csv.DictReader(open(args.config), delimiter='\t'):
        if rec['genome'] in args.targetGenomes:
            file_map[rec['genome']] = (rec['consensusGenePred'], rec['cgpGenePred'], rec['genomeFasta'],
                                       rec.get('cgpIntronBits', None))
    args.file_map = frozendict(file_map)
    return args


if __name__ == '__main__':
    args = parse_args()
    luigi.build([CgpConsensus(args)], local_scheduler=True, workers=args.localCores)
