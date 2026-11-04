#!/usr/bin/env python3
"""
Description: computes the clonal expansion/contraction statistics (Fisher's
exact test between every pair of timepoints, per subject/origin) shared by
template_overlap.qmd and template_discovery_brief.qmd. Both notebooks used to
run this computation independently and in full; this script runs it once so
both can load the results instead.
"""

import argparse
import itertools
import os

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

CLONE_COLS = ['junction_aa', 'v_call', 'j_call']


def build_clonotypes_df(concat_csv, sample_table, subject_col, origin_col,
                         timepoint_col, timepoint_order_col, timepoint_order,
                         alias_col):
    meta = pd.read_csv(sample_table, sep=',')

    # alias_col is optional; default to sample name if missing.
    if alias_col not in meta.columns:
        meta[alias_col] = meta['sample']

    # timepoint_order_col isn't a samplesheet column - compute and inject it
    # here. timepoint_order (an ordered comma-separated list, e.g.
    # "Base,Week4,EOT") ranks timepoints by position in that list; any
    # timepoint not listed sorts after the listed ones. When timepoint_order
    # is empty, timepoints are ranked numerically if every value parses as a
    # number, otherwise alphabetically.
    if timepoint_order_col not in meta.columns:
        unique_timepoints = meta[timepoint_col].dropna().unique().tolist()
        if timepoint_order:
            order_list = [t.strip() for t in timepoint_order.split(',')]
            rank_map = {t: i for i, t in enumerate(order_list)}
            unlisted = sorted(t for t in unique_timepoints if t not in rank_map)
            for i, t in enumerate(unlisted):
                rank_map[t] = len(order_list) + i
        else:
            try:
                rank_map = {t: r for r, t in enumerate(sorted(unique_timepoints, key=lambda x: float(x)))}
            except (TypeError, ValueError):
                rank_map = {t: r for r, t in enumerate(sorted(unique_timepoints, key=str))}
        meta[timepoint_order_col] = meta[timepoint_col].map(rank_map)

    concat_df = pd.read_csv(concat_csv, sep='\t')
    concat_df = concat_df.merge(
        meta[['sample', origin_col, subject_col, alias_col, timepoint_col, timepoint_order_col]],
        on='sample', how='left'
    )

    clonotypes_df = concat_df.copy()
    clonotypes_df = clonotypes_df.dropna(subset=[timepoint_col, origin_col])
    clonotypes_df['timepoint_rank'] = clonotypes_df[timepoint_order_col]
    clonotypes_df = clonotypes_df.sort_values('timepoint_rank')

    sample_total_counts = clonotypes_df.groupby([subject_col, origin_col, timepoint_col])['duplicate_count'].sum().reset_index()
    sample_total_counts.rename(columns={'duplicate_count': 'total_counts'}, inplace=True)
    clonotypes_df = pd.merge(clonotypes_df, sample_total_counts, on=[subject_col, origin_col, timepoint_col])
    clonotypes_df['frequency'] = clonotypes_df['duplicate_count'] / clonotypes_df['total_counts']

    return clonotypes_df


def run_fisher_vectorized(count_post, count_pre, total_post, total_pre, alt_hyp):
    # One-sided Fisher's exact test reduces to a single vectorized
    # hypergeometric call since total_pre/total_post are fixed margins.
    M = total_pre + total_post
    n = count_post + count_pre
    if alt_hyp == 'greater':
        return hypergeom.sf(count_post - 1, M, n, total_post)
    return hypergeom.cdf(count_post, M, n, total_post)


def compute_dynamics(clonotypes_df, subject_col, origin_col, timepoint_col,
                      direction, fc_thold, signif_thold):
    """direction: 'up' (expansion) or 'down' (contraction)."""
    merge_cols = CLONE_COLS + [subject_col, origin_col]

    all_significant_results = []
    comparison_details = []
    counts_list = []
    count_key = 'expanded_count' if direction == 'up' else 'contracted_count'

    for (subject, origin), subject_df in clonotypes_df.groupby([subject_col, origin_col]):
        timepoints = subject_df.sort_values('timepoint_rank')[timepoint_col].unique().tolist()
        if len(timepoints) < 2:
            continue

        freq_pivot = subject_df.pivot_table(index=CLONE_COLS, columns=timepoint_col, values='frequency', fill_value=0)
        counts_pivot = subject_df.pivot_table(index=CLONE_COLS, columns=timepoint_col, values='duplicate_count', fill_value=0)
        subject_total_counts = subject_df.groupby(timepoint_col)['duplicate_count'].sum()

        for t_pre, t_post in itertools.combinations(timepoints, 2):
            comp_df = pd.DataFrame({
                'count_pre': counts_pivot[t_pre],
                'count_post': counts_pivot[t_post],
                'freq_pre': freq_pivot[t_pre],
                'freq_post': freq_pivot[t_post]
            }).reset_index()

            if direction == 'up':
                # Must have post counts > 0 to be an expansion.
                comp_df = comp_df[comp_df['count_post'] > 0].copy()
            else:
                # Must have pre counts > 0 to contract.
                comp_df = comp_df[comp_df['count_pre'] > 0].copy()
            if comp_df.empty:
                continue

            comp_df['total_pre'] = subject_total_counts[t_pre]
            comp_df['total_post'] = subject_total_counts[t_post]

            if direction == 'up':
                comp_df['fold_change'] = np.where(comp_df['freq_pre'] == 0, np.inf, comp_df['freq_post'] / comp_df['freq_pre'])
                # Filter strictly by fold_change threshold before running the
                # expensive Fisher tests - if it doesn't meet the FC
                # threshold, it won't make the final cut anyway.
                comp_df = comp_df[comp_df['fold_change'] > fc_thold].copy()
            else:
                comp_df['fold_change'] = comp_df['freq_post'] / comp_df['freq_pre']
                comp_df = comp_df[comp_df['fold_change'] < fc_thold].copy()

            if comp_df.empty:
                continue

            comp_df['p_value'] = run_fisher_vectorized(
                comp_df['count_post'], comp_df['count_pre'], comp_df['total_post'], comp_df['total_pre'],
                'greater' if direction == 'up' else 'less'
            )

            comp_df[subject_col] = subject
            comp_df[origin_col] = origin
            comp_df['t_pre'] = t_pre
            comp_df['t_post'] = t_post

            _, q_values, _, _ = multipletests(comp_df['p_value'], alpha=signif_thold, method='fdr_bh')
            comp_df['q_value'] = q_values

            comparison_details.append(comp_df)

            significant_hits = comp_df[comp_df['q_value'] < signif_thold]

            if not significant_hits.empty:
                counts_per_pair = pd.DataFrame({
                    't_pre': [t_pre], 't_post': [t_post],
                    count_key: [len(significant_hits)],
                    subject_col: [subject], origin_col: [origin]
                })
                counts_list.append(counts_per_pair)

                significant_clones = significant_hits[CLONE_COLS].drop_duplicates()
                final_clones_df = freq_pivot.loc[pd.MultiIndex.from_frame(significant_clones)].reset_index()
                final_clones_df[subject_col] = subject
                final_clones_df[origin_col] = origin

                cols = [subject_col, origin_col] + CLONE_COLS + timepoints
                final_clones_df = final_clones_df[[c for c in cols if c in final_clones_df.columns]]
                all_significant_results.append(final_clones_df)

    # Empty results still carry their header: a patient with one timepoint per origin
    # (common - cohorts rarely share a sampling schedule) yields no comparisons, and a
    # header-less CSV is unreadable downstream.
    counts_table = pd.DataFrame(columns=[subject_col, origin_col, 't_pre', 't_post', count_key])
    if counts_list:
        counts_table = pd.concat(counts_list, ignore_index=True)
        counts_table = counts_table[[subject_col, origin_col, 't_pre', 't_post', count_key]]

    detailed_table = pd.DataFrame(columns=CLONE_COLS + [
        'count_pre', 'count_post', 'freq_pre', 'freq_post', 'total_pre', 'total_post',
        'fold_change', 'p_value', subject_col, origin_col, 't_pre', 't_post', 'q_value'])
    if comparison_details:
        detailed_table = pd.concat(comparison_details, ignore_index=True)
        # fold_change is only ever inf for expansions (freq_pre == 0, a brand
        # new clone) - contractions require count_pre > 0 so freq_pre can't be
        # 0. Both notebooks that consume this table treat these as "no finite
        # fold change to report" (NaN), not literal infinity.
        inf_cols = ['p_value', 'q_value', 'fold_change'] if direction == 'up' else ['p_value', 'q_value']
        for col in inf_cols:
            detailed_table[col] = detailed_table[col].replace([np.inf, -np.inf], np.nan)

    summary_table = pd.DataFrame(columns=[subject_col, origin_col] + CLONE_COLS)
    if all_significant_results:
        summary_table = pd.concat(all_significant_results, ignore_index=True).fillna(0)
        summary_table.rename(columns={t: f'freq_{t}' for t in summary_table.columns if t not in merge_cols}, inplace=True)

    return counts_table, detailed_table, summary_table


def main():
    parser = argparse.ArgumentParser(description="Compute clonal expansion/contraction dynamics")
    parser.add_argument('--concat_cdr3', required=True, help='concatenated_cdr3_sorted.tsv')
    parser.add_argument('--sample_table', required=True, help='samplesheet CSV')
    parser.add_argument('--subject_col', required=True)
    parser.add_argument('--origin_col', default='origin')
    parser.add_argument('--timepoint_col', required=True)
    parser.add_argument('--timepoint_order_col', required=True)
    parser.add_argument('--timepoint_order', default='')
    parser.add_argument('--alias_col', required=True)
    parser.add_argument('--fc_thold_up', type=float, default=2.0)
    parser.add_argument('--fc_thold_down', type=float, default=0.5)
    parser.add_argument('--signif_thold', type=float, default=0.05)
    parser.add_argument('--outdir', default='.')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    clonotypes_df = build_clonotypes_df(
        args.concat_cdr3, args.sample_table, args.subject_col, args.origin_col,
        args.timepoint_col, args.timepoint_order_col, args.timepoint_order, args.alias_col
    )
    clonotypes_df.to_csv(os.path.join(args.outdir, 'clonotypes.csv'), index=False)

    expansion_counts, detailed_expansion, expanded_summary = compute_dynamics(
        clonotypes_df, args.subject_col, args.origin_col, args.timepoint_col,
        'up', args.fc_thold_up, args.signif_thold
    )
    contraction_counts, detailed_contraction, contracting_summary = compute_dynamics(
        clonotypes_df, args.subject_col, args.origin_col, args.timepoint_col,
        'down', args.fc_thold_down, args.signif_thold
    )

    expansion_counts.to_csv(os.path.join(args.outdir, 'expansion_counts.csv'), index=False)
    detailed_expansion.to_csv(os.path.join(args.outdir, 'detailed_comparisons_expansion.csv'), index=False)
    expanded_summary.to_csv(os.path.join(args.outdir, 'expanded_clones_summary.csv'), index=False)

    contraction_counts.to_csv(os.path.join(args.outdir, 'contraction_counts.csv'), index=False)
    detailed_contraction.to_csv(os.path.join(args.outdir, 'detailed_comparisons_contraction.csv'), index=False)
    contracting_summary.to_csv(os.path.join(args.outdir, 'contracting_clones_summary.csv'), index=False)


if __name__ == "__main__":
    main()
