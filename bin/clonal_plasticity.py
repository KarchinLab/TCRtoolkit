#!/usr/bin/env python3
"""
CLONAL_PLASTICITY

Clone-level plasticity: how far a clone's cell-state distribution moves between timepoints.

For clone c and successive timepoints t1, t2, with state distributions P and Q over the cells
of that clone:

    plasticity(c, t1, t2) = JSD(P, Q),   base 2, so 0 <= JSD <= 1

Uncertainty: a percentile bootstrap resampling cells within each (clone, timepoint).

Small clones have positive JSD by chance alone, so the raw value is never read directly. Two
permutation nulls answer different questions:

time   shuffles timepoint labels among the clone's own cells. Asks whether this clone's state
       distribution changed over time. A population-wide shift also rejects this null.
state  shuffles state labels among all cells within each timepoint, preserving each
       timepoint's composition. Asks whether the clone moved more than a random set of cells
       of the same size drawn from the same timepoints, so it discounts population-wide shifts.

Inputs:
    --cells        per-cell export (sample, patient, clone_id, annot)
    --samplesheet  report samplesheet (sample, timepoint, origin)

Output (--outdir): clonal_plasticity.tsv, one row per patient, compartment, clone and pair of
successive timepoints, with q-values within each null.
"""
import argparse
import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def jsd(p, q):
    """Jensen-Shannon divergence in bits between count or probability vectors."""
    p = np.asarray(p, float)
    q = np.asarray(q, float)
    p = p / p.sum(axis=-1, keepdims=True)
    q = q / q.sum(axis=-1, keepdims=True)
    m = 0.5 * (p + q)

    def kl(a):
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(a > 0, a * np.log2(a / m), 0.0)
        return t.sum(axis=-1)

    return 0.5 * kl(p) + 0.5 * kl(q)


def _counts(codes, n_states):
    return np.bincount(codes, minlength=n_states)


def _rows_to_counts(matrix, n_states):
    """Per-row state counts for an (R, n) array of integer state codes."""
    offset = matrix + n_states * np.arange(matrix.shape[0])[:, None]
    return np.bincount(offset.ravel(), minlength=n_states * matrix.shape[0]).reshape(-1, n_states)


def clone_plasticity(cells, clone_col="clone", time_col="timepoint", state_col="state",
                     order=None, min_cells=5, n_boot=1000, n_perm=1000, seed=0):
    """Plasticity for each clone and pair of successive timepoints it was sampled at.

    Timepoints are visited in `order`, or sorted order. A pair qualifies when the clone has at
    least min_cells at both. Returns one row per (clone, t1, t2) with the JSD, a bootstrap
    95% interval, the mean of each null, and empirical p-values with the +1 correction.
    """
    rng = np.random.default_rng(seed)
    d = cells[[clone_col, time_col, state_col]].dropna()
    states, codes = np.unique(d[state_col].astype(str), return_inverse=True)
    d = d.assign(_s=codes)
    n_states = len(states)
    times = list(order) if order is not None else sorted(d[time_col].unique())
    state_totals = {t: _counts(g["_s"].to_numpy(), n_states) for t, g in d.groupby(time_col)}

    rows = []
    for clone, g in d.groupby(clone_col):
        cells_at = {t: x["_s"].to_numpy() for t, x in g.groupby(time_col)}
        ok = [t for t in times if len(cells_at.get(t, ())) >= min_cells]
        for t1, t2 in zip(ok, ok[1:]):
            a, b = cells_at[t1], cells_at[t2]
            obs = float(jsd(_counts(a, n_states), _counts(b, n_states)))

            boot = jsd(_rows_to_counts(rng.choice(a, (n_boot, len(a))), n_states),
                       _rows_to_counts(rng.choice(b, (n_boot, len(b))), n_states))

            pooled = np.concatenate([a, b])
            perm = rng.permuted(np.tile(pooled, (n_perm, 1)), axis=1)
            null_time = jsd(_rows_to_counts(perm[:, :len(a)], n_states),
                            _rows_to_counts(perm[:, len(a):], n_states))

            # Exact draws of len(a) and len(b) cells without replacement from each timepoint.
            null_state = jsd(rng.multivariate_hypergeometric(state_totals[t1], len(a), size=n_perm),
                             rng.multivariate_hypergeometric(state_totals[t2], len(b), size=n_perm))

            rows.append({"clone": clone, "t1": t1, "t2": t2, "n1": len(a), "n2": len(b),
                         "jsd": obs, "ci_low": float(np.quantile(boot, 0.025)),
                         "ci_high": float(np.quantile(boot, 0.975)),
                         "null_time_mean": float(null_time.mean()),
                         "null_state_mean": float(null_state.mean()),
                         "p_time": float((1 + (null_time >= obs - 1e-12).sum()) / (n_perm + 1)),
                         "p_state": float((1 + (null_state >= obs - 1e-12).sum()) / (n_perm + 1))})
    return pd.DataFrame(rows)


RESULT_COLS = ["patient", "origin", "clone", "t1", "t2", "n1", "n2", "jsd", "ci_low", "ci_high",
               "null_time_mean", "null_state_mean", "p_time", "p_state", "q_time", "q_state"]


def timepoint_order(values, explicit):
    """Explicit order, else numeric when every value parses, else alphabetical."""
    values = sorted(set(values))
    if explicit:
        wanted = [x.strip() for x in explicit.split(",") if x.strip()]
        return [t for t in wanted if t in values] + [t for t in values if t not in wanted]
    numeric = pd.to_numeric(pd.Series(values), errors="coerce")
    return [values[i] for i in numeric.argsort()] if numeric.notna().all() else values


def run(cells, sheet, order, min_cells, n_boot, n_perm, seed, outdir):
    os.makedirs(outdir, exist_ok=True)
    meta = sheet.drop_duplicates("sample").set_index("sample")
    d = cells[["sample", "patient", "clone_id", "annot"]].copy()
    d["timepoint"] = d["sample"].map(meta["timepoint"]) if "timepoint" in meta else np.nan
    d["origin"] = d["sample"].map(meta["origin"]).fillna("all") if "origin" in meta else "all"
    d = d.dropna(subset=["patient", "clone_id", "annot", "timepoint"])
    d["timepoint"] = d["timepoint"].astype(str)
    times = timepoint_order(d["timepoint"], order)
    parts = []
    for (patient, origin), g in d.groupby(["patient", "origin"]):
        r = clone_plasticity(g, "clone_id", "timepoint", "annot", times, min_cells, n_boot, n_perm, seed)
        if len(r):
            parts.append(r.assign(patient=patient, origin=origin))
    out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=RESULT_COLS)
    if len(out):
        for col in ("time", "state"):
            out[f"q_{col}"] = multipletests(out[f"p_{col}"], method="fdr_bh")[1]
    out[RESULT_COLS].to_csv(os.path.join(outdir, "clonal_plasticity.tsv"), sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser(description="Clone-level state change between timepoints")
    ap.add_argument("--cells", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--timepoint-order", default="")
    ap.add_argument("--min-cells", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    run(pd.read_csv(a.cells, sep="\t", low_memory=False), pd.read_csv(a.samplesheet), a.timepoint_order,
        a.min_cells, a.n_boot, a.n_perm, a.seed, a.outdir)


if __name__ == "__main__":
    main()
