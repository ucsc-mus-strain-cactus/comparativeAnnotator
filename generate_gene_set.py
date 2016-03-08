"""
Produces a gene set from transMap alignments, or from a combination of transMap and AugustusTM/TMR.
"""
import cPickle as pickle
from collections import defaultdict, OrderedDict
from comparativeAnnotator.database_queries import get_row_dict, get_fail_pass_excel_ids, augustus_eval
from pycbio.sys.dataOps import merge_dicts
from pycbio.sys.mathOps import format_ratio
from pycbio.bio.transcripts import GenePredTranscript, get_transcript_dict
from pycbio.bio.intervals import ChromosomeInterval
from comparativeAnnotator.comp_lib.name_conversions import strip_alignment_numbers, remove_augustus_alignment_number, \
    aln_id_is_augustus, aln_id_is_transmap
from comparativeAnnotator.database_queries import get_gene_transcript_map, get_transcript_gene_map, \
    get_transcript_biotype_map

__author__ = "Ian Fiddes"


def mode_is_aug(mode):
    assert mode in ['AugustusTMR', 'AugustusTM', 'transMap', 'augustus']
    return mode == 'AugustusTMR' or mode == 'AugustusTM' or mode == 'augustus'


def load_gps(gp_paths):
    """
    Get a dictionary mapping all gene IDs from a genePred into its entire record. If the gene IDs are not unique
    this function will not work like you want it to.
    """
    r = {}
    for gp in gp_paths:
        r.update(get_transcript_dict((gp)))
    return r


def get_db_rows(ref_genome, genome, db, biotype, mode):
    """
    Adapter function to return the combination of the database rows produced for augustus and transMap
    """
    tm_stats = get_row_dict(ref_genome, genome, db, 'transMap', biotype)
    if mode_is_aug(mode):
        aug_stats = get_row_dict(ref_genome, genome, db, mode, biotype)
        return merge_dicts([tm_stats, aug_stats])
    else:
        return tm_stats


def build_data_dict(id_names, id_list, transcript_gene_map, gene_transcript_map):
    """
    Builds a dictionary mapping gene_id -> transcript_ids -> aln_ids in id_names bins (as an OrderedDict)
    """
    data_dict = defaultdict(dict)
    for gene_id in gene_transcript_map:
        for ens_id in gene_transcript_map[gene_id]:
            data_dict[gene_id][ens_id] = OrderedDict((x, []) for x in id_names)
    for ids, n in zip(*[id_list, id_names]):
        for aln_id in ids:
            ens_id = strip_alignment_numbers(aln_id)
            gene_id = transcript_gene_map[ens_id]
            if gene_id in data_dict and ens_id in data_dict[gene_id]:
                data_dict[gene_id][ens_id][n].append(aln_id)
    return data_dict


def find_best_alns(stats, tm_ids, aug_ids, tm_cov_cutoff, aug_cov_cutoff):
    """
    Takes the list of transcript Ids and finds the best alignment(s) by highest percent identity and coverage
    We sort by ID to favor Augustus transcripts going to the consensus set in the case of ties
    """
    def get_cov_ident(stats, aln_id, mode):
        """extract coverage and identity from stats, round them"""
        if mode == 'transMap':
            cov = stats[aln_id].AlignmentCoverage
            ident = stats[aln_id].AlignmentIdentity
        elif mode_is_aug(mode):
            cov = stats[aln_id].AugustusAlignmentCoverage
            ident = stats[aln_id].AugustusAlignmentIdentity
        else:
            raise NotImplementedError
        # round to avoid floating point issues when finding ties
        if cov is None:
            cov = 0.0
        else:
            cov = round(cov, 6)
        if ident is None:
            ident = 0.0
        else:
            ident = round(ident, 6)
        return cov, ident
    def analyze_ids(ids, cutoff, stats, mode):
        """return all ids which pass coverage cutoff"""
        s = []
        for aln_id in ids:
            cov, ident = get_cov_ident(stats, aln_id, mode)
            s.append([aln_id, cov, ident])
        return filter(lambda (aln_id, cov, ident): cov >= cutoff, s)
    aug_cov = analyze_ids(aug_ids, aug_cov_cutoff, stats, mode='augustus')
    tm_cov = analyze_ids(tm_ids, tm_cov_cutoff, stats, mode='transMap')
    cov_s = aug_cov + tm_cov
    if len(cov_s) == 0:
        return None
    else:
        best_ident = sorted(cov_s, key=lambda (aln_id, cov, ident): ident, reverse=True)[0][2]
    best_overall = [aln_id for aln_id, cov, ident in cov_s if ident >= best_ident]
    return best_overall


def evaluate_ids(fail_ids, pass_specific_ids, excel_ids, aug_ids, stats, tm_cov_cutoff, aug_cov_cutoff):
    """
    For a given ensembl ID, we have augustus/transMap ids in 3 categories. Based on the hierarchy Excellent>Pass>Fail,
    return the best transcript in the highest category with a transMap transcript.
    """
    if len(excel_ids) > 0:
        best_alns = find_best_alns(stats, excel_ids, aug_ids, tm_cov_cutoff, aug_cov_cutoff)
        return best_alns, "Excellent"
    elif len(pass_specific_ids) > 0:
        best_alns = find_best_alns(stats, pass_specific_ids, aug_ids, tm_cov_cutoff, aug_cov_cutoff)
        return best_alns, "Pass"
    elif len(fail_ids) > 0:
        best_alns = find_best_alns(stats, fail_ids, aug_ids, tm_cov_cutoff, aug_cov_cutoff)
        return best_alns, "Fail"
    else:
        return None, "NoTransMap"


def is_tie(best_alns):
    """
    If we have more than one best transcript, is at least one from transMap and one from Augustus?
    """
    seen = set()
    for aln_id in best_alns:
        ens_id = remove_augustus_alignment_number(aln_id)
        if ens_id in seen:
            return True
        else:
            seen.add(ens_id)
    return False


def find_best_transcripts(data_dict, stats, mode, biotype, tm_cov_cutoff=80.0, aug_cov_cutoff=50.0):
    """
    For all of the transcripts categorized in data_dict, evaluate them and bin them.
    """
    binned_transcripts = {}
    for gene_id in data_dict:
        binned_transcripts[gene_id] = {}
        for ens_id in data_dict[gene_id]:
            tx_recs = data_dict[gene_id][ens_id]
            if mode_is_aug(mode) and biotype == "protein_coding":
                fail_ids, pass_specific_ids, excel_ids, aug_ids = tx_recs.values()
            else:
                fail_ids, pass_specific_ids, excel_ids = tx_recs.values()
                aug_ids = []
            best_alns, category = evaluate_ids(fail_ids, pass_specific_ids, excel_ids, aug_ids, stats, tm_cov_cutoff,
                                               aug_cov_cutoff)
            if best_alns is None:
                binned_transcripts[gene_id][ens_id] = [best_alns, category, None]
            else:
                tie = is_tie(best_alns)
                binned_transcripts[gene_id][ens_id] = [best_alns[0], category, tie]
    return binned_transcripts


def find_longest_for_gene(bins, stats, gps, cov_cutoff=33.3, ident_cutoff=80.0):
    """
    Finds the longest transcript(s) for a gene. This is used when all transcripts failed, and has more relaxed cutoffs.
    """
    def tm_has_gaps(tm_stats):
        return tm_stats.UtrGap + tm_stats.CdsGap + tm_stats.UnknownGap + tm_stats.CdsMult3Gap != 0
    aln_ids = zip(*bins.itervalues())[0]
    keep_ids = []
    for aln_id in aln_ids:
        if aln_id is None:
            continue
        elif aln_id_is_transmap(aln_id):
            # we don't allow transMap transcripts through with short introns
            if tm_has_gaps(stats[aln_id]):
                continue
            cov = stats[aln_id].AlignmentCoverage
            ident = stats[aln_id].AlignmentIdentity
        elif aln_id_is_augustus(aln_id):
            cov = stats[aln_id].AugustusAlignmentCoverage
            ident = stats[aln_id].AugustusAlignmentIdentity
        else:
            raise NotImplementedError
        if cov >= cov_cutoff and ident >= ident_cutoff:
            keep_ids.append(aln_id)
    if len(keep_ids) > 0:
        sizes = [[x, len(gps[x])] for x in keep_ids]
        longest_size = max(zip(*sizes)[1])
        r = [x for x, y in sizes if y == longest_size]
        return r[0], is_tie(r)
    else:
        return None, None


def has_only_short(ids_included, ref_size, gps, percentage_of_ref=70.0):
    """
    Are all of the consensus transcripts we found for this gene too short?
    """
    r = []
    for aln_id in ids_included:
        tgt_size = len(gps[aln_id])
        r.append(tgt_size)
    return all([100 * format_ratio(tgt_size, ref_size) < percentage_of_ref for tgt_size in r])


def evaluate_consensus_tx(best_id, category, tie):
    """
    Evaluates the best transcript(s) for a given ensembl ID for being excel/fail/ok and asks if it is a tie
    """
    if tie is True:
        c = "Tie"
    elif aln_id_is_augustus(best_id):
        c = "Aug"
    elif aln_id_is_transmap(best_id):
        c = "TM"
    else:
        assert False, "ID was not TM/Aug"
    s = "".join([category, c])
    return s


def evaluate_gene(categories):
    """
    Same as evaluate_transcript, but on the gene level. Does this gene have at least one transcript categorized
    as excellent/passing/fail?
    """
    if "Excellent" in categories:
        return "Excellent"
    elif "Pass" in categories:
        return "Pass"
    elif "Fail" in categories:
        return "Fail"
    elif "NoTransMap" in categories:
        return "NoTransMap"
    else:
        assert False, "Should not be able to get here."


def generate_gene_set(binned_transcripts, stats, gps, transcript_biotype_map, ref_gene_sizes, mode, biotype):
    """
    Takes the binned transcripts and builds a consensus gene set.
    """
    if mode_is_aug(mode) and biotype == "protein_coding":
        is_consensus = True
        transcript_evaluation = OrderedDict((x, 0) for x in ["ExcellentTM", "ExcellentAug", "ExcellentTie",
                                                             "PassTM", "PassAug", "PassTie",
                                                             "Fail"])
    else:
        is_consensus = False
        transcript_evaluation = OrderedDict((x, 0) for x in ["Excellent", "Pass", "Fail"])
    gene_evaluation = OrderedDict((x, 0) for x in ["Excellent", "Pass", "Fail", "NoTransMap"])
    longest_rate = OrderedDict((('Longest', OrderedDict((('AddLongest', 0), ('FailAddLongest', 0)))),
                                ('Rescue', OrderedDict((('GeneRescue', 0), ('FailGeneRescue', 0))))))
    consensus = []
    for gene_id in binned_transcripts:
        ids_included = set()
        categories = set()
        for ens_id in binned_transcripts[gene_id]:
            # evaluate each transcript for a gene
            best_id, category, tie = binned_transcripts[gene_id][ens_id]
            categories.add(category)
            if category in ["Excellent", "Pass"]:
                # if a transcript is of high quality, include it
                consensus.append(best_id)
                ids_included.add(best_id)
                s = evaluate_consensus_tx(best_id, category, tie) if is_consensus is True else category
                transcript_evaluation[s] += 1
        # have we included this gene yet? only count if the transcripts included match the gene biotype.
        # this prevents good mappings of retained introns and such being the only transcript.
        biotype_ids_included = {x for x in ids_included if transcript_biotype_map[strip_alignment_numbers(x)] == biotype}
        gene_in_consensus = True if len(biotype_ids_included) > 0 else False
        # if we have, have we included only short transcripts?
        has_only_short_txs = has_only_short(ids_included, ref_gene_sizes[gene_id], gps)
        if has_only_short_txs is True and gene_in_consensus is True:
            # add the single longest transcript for this gene if it passes filters
            longest_id, tie = find_longest_for_gene(binned_transcripts[gene_id], stats, gps)
            if longest_id is not None and longest_id not in consensus:
                consensus.append(longest_id)
                transcript_evaluation['Fail'] += 1
                longest_rate['Longest']["AddLongest"] += 1
            else:
                longest_rate['Longest']["FailAddLongest"] += 1
        if gene_in_consensus is True:
            # gene is in consensus, evaluate and move on
            s = evaluate_gene(categories)
            gene_evaluation[s] += 1
        else:
            # attempt to add one longest transcript for this failing gene
            longest_id, tie = find_longest_for_gene(binned_transcripts[gene_id], stats, gps)
            if longest_id is None:
                gene_evaluation['NoTransMap'] += 1
                longest_rate['Rescue']["FailGeneRescue"] += 1
            else:
                category = 'Fail'
                consensus.append(longest_id)
                transcript_evaluation[category] += 1
                gene_evaluation[category] += 1
                longest_rate['Rescue']["GeneRescue"] += 1
    metrics = {"transcript": transcript_evaluation, "gene": gene_evaluation, "longest": longest_rate}
    return consensus, metrics


def consensus_by_biotype(db_path, ref_genome, genome, biotype, gps, transcript_gene_map, gene_transcript_map,
                         transcript_biotype_map, stats, mode, ref_gene_sizes, filter_chroms):
    """
    Main consensus finding function.
    """
    excel_ids, pass_specific_ids, fail_ids = get_fail_pass_excel_ids(ref_genome, genome, db_path, biotype,
                                                                     filter_chroms, best_cov_only=False)
    # hacky way to avoid duplicating code in consensus finding - we will always have an aug_id set, it just may be empty
    if mode_is_aug(mode) and biotype == "protein_coding":
        aug_ids = augustus_eval(ref_genome, genome, db_path, biotype, filter_chroms)
        id_names = ["fail_ids", "pass_specific_ids", "excel_ids", "aug_ids"]
        id_list = [fail_ids, pass_specific_ids, excel_ids, aug_ids]
    else:
        id_names = ["fail_ids", "pass_specific_ids", "excel_ids"]
        id_list = [fail_ids, pass_specific_ids, excel_ids]
    data_dict = build_data_dict(id_names, id_list, transcript_gene_map, gene_transcript_map)
    binned_transcripts = find_best_transcripts(data_dict, stats, mode, biotype)
    consensus, metrics = generate_gene_set(binned_transcripts, stats, gps, transcript_biotype_map, ref_gene_sizes, 
                                           mode, biotype)
    return consensus, metrics


def deduplicate_consensus(consensus, gps, stats):
    """
    In the process of consensus building, we may find that we have ended up with more than one transcript for a gene
    that are actually identical. Remove these, picking the best based on the stats dict.
    """
    duplicates = defaultdict(list)
    for tx_id in consensus:
        tx = gps[tx_id]
        duplicates[frozenset(tx.exon_intervals)].append(tx)
    deduplicated_consensus = []
    dup_count = 0
    for gp_list in duplicates.itervalues():
        if len(gp_list) > 1:
            dup_count += 1
            # we have duplicates to collapse - which has the highest %ID followed by highest %coverage?
            dup_stats = sorted([[x, stats[x.name]] for x in gp_list], key=lambda (n, r): (r.AlignmentIdentity,
                                                                                          r.AlignmentCoverage))
            best = dup_stats[0][0].name
            deduplicated_consensus.append(best)
        else:
            deduplicated_consensus.append(gp_list[0].name)
    return deduplicated_consensus, dup_count


def fix_gene_pred(gp, transcript_gene_map):
    """
    These genePreds have a few problems. First, the alignment numbers must be removed. Second, we want to fix
    the name2 field to be the gene name. Third, we want to set the unique ID field. Finally, we want to sort the whole
    thing by genomic coordinates.
    Also reports the number of genes and transcripts seen.
    """
    gp = sorted(gp, key=lambda tx: (tx.chromosome, tx.start))
    fixed = []
    for tx in gp:
        x = tx.get_gene_pred()
        x[10] = x[0]  # use unique Aug/TM ID as unique identifier
        tx_id = strip_alignment_numbers(x[0])
        x[0] = tx_id
        gene_id = transcript_gene_map[tx_id]
        x[11] = gene_id
        fixed.append(x)
    return ["\t".join(map(str, x)) for x in fixed]


def write_gps(consensus, gps, gp_path, transcript_gene_map):
    """
    Writes the final consensus gene set to a genePred, after fixing the names. Reports the number of genes and txs
    in the final set
    """
    gp_recs = [gps[aln_id] for aln_id in consensus]
    fixed_gp_recs = fix_gene_pred(gp_recs, transcript_gene_map)
    with open(gp_path, "w") as outf:
        for rec in fixed_gp_recs:
            outf.write(rec + '\n')


def build_gene_sizes(tx_dict, gene_transcript_map):
    """
    Finds the largest size transcript for a gene
    """
    r = {}
    for gene_id, tx_ids in gene_transcript_map.iteritems():
        r[gene_id] = max([len(tx_dict[x]) for x in tx_ids])
    return r


def generate_gene_set_wrapper(args):
    assert args.mode in ['AugustusTMR', 'AugustusTM', 'transMap']
    if mode_is_aug(args.mode):
        gps = load_gps([args.gp, args.augustus_gp])
    else:
        gps = load_gps([args.gp])
    transcript_gene_map = get_transcript_gene_map(args.query_genome, args.db)
    transcript_biotype_map = get_transcript_biotype_map(args.query_genome, args.db)
    ref_gps = get_transcript_dict(args.annotation_gp)
    biotype_evals = {}
    for biotype, gp_path in args.geneset_gps.iteritems():
        gene_transcript_map = get_gene_transcript_map(args.query_genome, args.db, biotype)
        ref_gene_sizes = build_gene_sizes(ref_gps, gene_transcript_map)
        stats = get_db_rows(args.query_genome, args.target_genome, args.db, biotype, args.mode)
        consensus, metrics = consensus_by_biotype(args.db, args.query_genome, args.target_genome, biotype, gps,
                                                  transcript_gene_map, gene_transcript_map, transcript_biotype_map,
                                                  stats, args.mode, ref_gene_sizes, args.filter_chroms)
        deduplicated_consensus, dup_count = deduplicate_consensus(consensus, gps, stats)
        write_gps(consensus, gps, gp_path, transcript_gene_map)
        metrics["duplication_rate"] = dup_count
        biotype_evals[biotype] = metrics
    with open(args.pickled_metrics, 'w') as outf:
        pickle.dump(biotype_evals, outf)
