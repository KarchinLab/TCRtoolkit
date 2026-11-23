#!/usr/bin/env python3
"""
SELECTION_TEST

Generative-null test for clonotype sharing.

Model
-----
Subject s contributes N_s distinct CDR3beta amino-acid clonotypes. Under recombination only,
clonotype i with OLGA generation probability p_i is present in subject s with probability

    pi_is(q) = 1 - exp(-c(p_i) * q * p_i * N_s)

c(p)  Calibration curve. OLGA's P_gen is not normalised to the productive sequences a
      repertoire contains, and c(p) also absorbs the mean effect of thymic selection at
      that P_gen. It is estimated from an anchor sample drawn from the generative model
      itself, whose cohort sharing K (including K = 0) is unbiased by observation.
q     Per-sequence heterogeneity at fixed P_gen, Gamma distributed with mean 1 and variance
      phi(p). phi(p) = 0 recovers the pure recombination null.

K_i, the number of subjects sharing clonotype i, is Poisson-binomial given q. Only observed
clonotypes are tested, so the p-value conditions on K_i >= 1:

    P_i = P(K >= k_i | K >= 1)

Tail probabilities are tabulated once on a grid of log10 rates and shared by all sequences.

Inputs:
    --clones      concatenated clone table with P_gen (junction_aa, sample, pgen)
    --samplesheet sample sheet; the subject column defines who counts as one subject
    --anchors     assets/olga/anchors_human_trb.tsv.gz (cdr3, pgen): sequences drawn from the
                  OLGA model, made once by assets/olga/make_anchors.py

Outputs (--outdir):
    selection_test.tsv         clonotypes in >= 2 subjects: observed and expected sharing,
                               P_gen, p, q and the call
    selection_summary.tsv      metric / value
    selection_calibration.tsv  the fitted c(p) and phi(p)
"""
import argparse
import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import optimize, stats

MAX_LOG10_PHI = 1.5


def poisson_binomial_pmf(probs):
    """Exact PMF of a sum of independent Bernoullis, for many rows of probabilities.

    probs has shape (R, S); returns (R, S + 1).
    """
    probs = np.atleast_2d(probs)
    pmf = np.zeros((probs.shape[0], probs.shape[1] + 1))
    pmf[:, 0] = 1.0
    for j in range(probs.shape[1]):
        p = probs[:, j:j + 1]
        pmf[:, 1:j + 2] = pmf[:, 1:j + 2] * (1 - p) + pmf[:, :j + 1] * p
        pmf[:, 0:1] *= (1 - p)
    return pmf


def presence_prob(rate, depths):
    """pi for each rate (R,) and subject depth (S,), as an (R, S) array."""
    return -np.expm1(-np.outer(np.asarray(rate, float), np.asarray(depths, float)))


def gamma_nodes(phi, n=24):
    """Equal-weight quantile nodes of Gamma(mean 1, variance phi), rescaled to mean 1."""
    if phi <= 1e-8:
        return np.ones(1)
    u = (np.arange(n) + 0.5) / n
    q = stats.gamma.ppf(u, a=1.0 / phi, scale=phi)
    return q / q.mean()


@dataclass
class Calibration:
    """Piecewise-linear c(p) and phi(p) on log10 P_gen."""
    log10_pgen: np.ndarray
    log10_c: np.ndarray
    phi: np.ndarray

    def c(self, pgen):
        return 10 ** np.interp(np.log10(pgen), self.log10_pgen, self.log10_c)

    def phi_at(self, pgen):
        return np.interp(np.log10(pgen), self.log10_pgen, self.phi)

    def to_frame(self):
        return pd.DataFrame({"log10_pgen": self.log10_pgen, "log10_c": self.log10_c, "phi": self.phi})


def calibrate(anchor_pgen, anchor_k, depths, n_bins=12, trim=0.999, fit_phi=True, min_events=30,
              max_moment_seqs=3000):
    """Fit c(p) and phi(p) from anchor sequences drawn from the generative model.

    anchor_k is the number of subjects containing each anchor, zeros included. Sharing is
    winsorised at the `trim` quantile within each bin so that the rare genuinely selected
    anchors do not set the null. Bins with fewer than `min_events` positive anchors are
    merged into their neighbour. Within a bin, c and phi jointly match the mean and variance
    of K with q integrated out; fitting them with q fixed at 1 biases both where presence
    saturates at high P_gen.
    """
    pgen = np.asarray(anchor_pgen, float)
    k = np.asarray(anchor_k, float)
    ok = pgen > 0
    pgen, k = pgen[ok], k[ok]
    depths = np.asarray(depths, float)
    order = np.argsort(pgen)
    pgen, k = pgen[order], k[order]

    # Equal-event bins: information about c comes from anchors seen at least once.
    cum = np.cumsum(k > 0)
    edges = np.searchsorted(cum, np.linspace(0, cum[-1], n_bins + 1)[1:-1])
    groups = [g for g in np.split(np.arange(len(k)), np.unique(edges)) if len(g)]
    merged = []
    for g in groups:
        if merged and (k[merged[-1]] > 0).sum() < min_events:
            merged[-1] = np.concatenate([merged[-1], g])
        else:
            merged.append(g)
    if len(merged) > 1 and (k[merged[-1]] > 0).sum() < min_events:
        merged[-2] = np.concatenate([merged[-2], merged[-1]])
        merged.pop()

    moments = _RateMoments(depths)

    rng = np.random.default_rng(0)
    xs, cs, phis = [], [], []
    for g in merged:
        p, kk = pgen[g], k[g]
        kk = np.minimum(kk, max(np.quantile(kk, trim), 1.0))
        target_mean, target_var = kk.mean(), kk.var()
        sub = p if len(p) <= max_moment_seqs else rng.choice(p, max_moment_seqs, replace=False)
        fit = lambda phi: _fit_c(sub, moments, phi, target_mean)
        phi = 0.0
        log_c, var0 = fit(0.0)
        if fit_phi and var0 < target_var:
            # Variance is not monotone in phi at extreme phi, so take the first crossing.
            ladder = np.arange(-3.0, MAX_LOG10_PHI + 0.01, 0.25)
            below = ladder[0]
            for x in ladder:
                if fit(10 ** x)[1] >= target_var:
                    break
                below = x
            lo, hi = below, x
            for _ in range(6):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if fit(10 ** mid)[1] < target_var else (lo, mid)
            phi = 10 ** (0.5 * (lo + hi))
            log_c = fit(phi)[0]
        xs.append(np.median(np.log10(p)))
        cs.append(log_c)
        phis.append(phi)
    return Calibration(np.array(xs), np.array(cs), np.array(phis))


class _RateMoments:
    """Mean and variance of K given a per-sequence rate, tabulated on a log10 rate grid."""

    def __init__(self, depths, lo=-22.0, hi=2.0, step=0.02):
        self.grid = np.arange(lo, hi + step, step)
        mean = np.zeros_like(self.grid)
        var = np.zeros_like(self.grid)
        for i in range(0, len(self.grid), 256):
            pi = presence_prob(10 ** self.grid[i:i + 256], depths)
            mean[i:i + 256] = pi.sum(axis=1)
            var[i:i + 256] = (pi * (1 - pi)).sum(axis=1)
        self.log_mean = np.log(np.maximum(mean, 1e-300))
        self.log_var = np.log(np.maximum(var, 1e-300))

    def __call__(self, log_rate):
        return (np.exp(np.interp(log_rate, self.grid, self.log_mean)),
                np.exp(np.interp(log_rate, self.grid, self.log_var)))

    def moments(self, log_p, nodes, log_c):
        """Mean and variance of K over sequences, with q integrated out."""
        log_rate = log_c + log_p[None, :] + np.log10(np.maximum(nodes, 1e-300))[:, None]
        mean_q, var_q = self(log_rate)
        mean = mean_q.mean()
        return mean, (var_q + mean_q ** 2).mean() - mean ** 2


def _fit_c(p, moments, phi, target_mean):
    """log10 c matching the mean sharing at this phi; returns (log10 c, implied variance)."""
    nodes = gamma_nodes(phi)
    log_p = np.log10(p)

    def excess(log_c):
        return moments.moments(log_p, nodes, log_c)[0] - target_mean

    if excess(-6) >= 0:
        log_c = -6.0
    elif excess(6) <= 0:
        log_c = 6.0
    else:
        log_c = optimize.brentq(excess, -6, 6, xtol=1e-3)
    return log_c, moments.moments(log_p, nodes, log_c)[1]


class SharingNull:
    """Tail table of P(K >= k) over a grid of log10 per-sequence rates."""

    def __init__(self, depths, grid_step=0.05, lo=-16.0, hi=-2.0):
        self.depths = np.asarray(depths, float)
        self.grid = np.arange(lo, hi + grid_step, grid_step)
        rates = 10 ** self.grid
        pmf = np.concatenate([
            poisson_binomial_pmf(presence_prob(rates[i:i + 64], self.depths))
            for i in range(0, len(rates), 64)])
        self.pmf = pmf
        self.tail = np.clip(np.cumsum(pmf[:, ::-1], axis=1)[:, ::-1], 0, 1)

    def _interp(self, table, log_rate, k):
        x = np.clip((log_rate - self.grid[0]) / (self.grid[1] - self.grid[0]), 0, len(self.grid) - 1.000001)
        i = x.astype(int)
        w = x - i
        lo = np.log(np.maximum(table[i, k], 1e-300))
        hi = np.log(np.maximum(table[i + 1, k], 1e-300))
        return np.exp(lo * (1 - w) + hi * w)

    def tail_prob(self, log_rate, k):
        return self._interp(self.tail, log_rate, k)

    def point_prob(self, log_rate, k):
        return self._interp(self.pmf, log_rate, k)


def expected_sharing(pgen, depths, calibration):
    rate = calibration.c(pgen) * pgen
    return presence_prob(rate, depths).sum(axis=1)


def sharing_pvalues(pgen, k, null, calibration, n_nodes=24, phi_bins=40):
    """Return (p, mid_p) for P(K >= k | K >= 1) under the calibrated null.

    p is the conservative tail, and mid_p counts half of the mass at k.
    Sequences with P_gen = 0 or missing return NaN: the model cannot generate them, so
    the null makes no statement about them.
    """
    pgen = np.asarray(pgen, float)
    k = np.asarray(k, int)
    p_out = np.full(len(pgen), np.nan)
    mid_out = np.full(len(pgen), np.nan)
    ok = np.isfinite(pgen) & (pgen > 0) & (k >= 1)
    if not ok.any():
        return p_out, mid_out
    pg, kk = pgen[ok], k[ok]
    log_rate = np.log10(calibration.c(pg) * pg)
    phi = calibration.phi_at(pg)

    # Sequences with similar phi share quadrature nodes.
    ladder = np.unique(np.quantile(phi, np.linspace(0, 1, phi_bins + 1)))
    bin_of = np.clip(np.searchsorted(ladder, phi, side="right") - 1, 0, max(len(ladder) - 2, 0))
    tail_k = np.zeros(len(pg))
    point_k = np.zeros(len(pg))
    tail_1 = np.zeros(len(pg))
    for b in np.unique(bin_of):
        m = bin_of == b
        nodes = gamma_nodes(float(phi[m].mean()), n_nodes)
        for q in nodes:
            lr = log_rate[m] + np.log10(max(q, 1e-300))
            tail_k[m] += null.tail_prob(lr, kk[m]) / len(nodes)
            point_k[m] += null.point_prob(lr, kk[m]) / len(nodes)
            tail_1[m] += null.tail_prob(lr, np.ones(m.sum(), int)) / len(nodes)
    denom = np.maximum(tail_1, 1e-300)
    p_out[ok] = np.minimum(tail_k / denom, 1.0)
    mid_out[ok] = np.minimum((tail_k - 0.5 * point_k) / denom, 1.0)
    return p_out, mid_out


def randomized_pit(pgen, k, null, calibration, rng, n_nodes=24):
    """Randomised probability integral transform, exactly uniform under a correct null.

    Discrete tests are never exactly uniform, so calibration is judged on this quantity.
    """
    p, mid = sharing_pvalues(pgen, k, null, calibration, n_nodes)
    point = 2 * (p - mid)
    return p - rng.random(len(p)) * point


def bh_adjust(p, m_total=None):
    """Benjamini-Hochberg q-values. m_total counts hypotheses absent from p, such as
    observed clonotypes that were never scored because they cannot reach significance."""
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    m = int(m_total) if m_total else int(ok.sum())
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx], kind="stable")]
    ranks = np.arange(1, len(order) + 1)
    adj = np.minimum.accumulate((p[order] * m / ranks)[::-1])[::-1]
    q[order] = np.minimum(adj, 1.0)
    return q


RESULT_COLS = ["junction_aa", "k", "expected_k", "pgen", "log10_pgen", "p_value", "q_value", "call"]


def subject_table(clones, samplesheet, subject_col):
    """(subject, junction_aa, pgen) rows. Falls back to one subject per sample, reported in
    the summary, because repeat samples of a person are not independent subjects."""
    sheet = pd.read_csv(samplesheet)
    unit = "sample"
    subject = clones["sample"]
    if subject_col in sheet.columns and sheet[subject_col].notna().any():
        mapped = clones["sample"].map(sheet.drop_duplicates("sample").set_index("sample")[subject_col])
        if mapped.notna().any():
            unit, subject = subject_col, mapped
    return clones.assign(subject=subject).dropna(subset=["subject", "junction_aa"]), unit


def run(clones, samplesheet, anchors, subject_col, fdr, min_subjects, outdir):
    os.makedirs(outdir, exist_ok=True)
    d, unit = subject_table(clones, samplesheet, subject_col)
    depths = d.groupby("subject")["junction_aa"].nunique()
    summary = [("Subject unit", unit), ("Subjects", len(depths)), ("Unique clonotypes", d["junction_aa"].nunique())]

    def finish(result, calibration):
        result.to_csv(os.path.join(outdir, "selection_test.tsv"), sep="\t", index=False)
        calibration.to_csv(os.path.join(outdir, "selection_calibration.tsv"), sep="\t", index=False)
        pd.DataFrame(summary, columns=["metric", "value"]).to_csv(
            os.path.join(outdir, "selection_summary.tsv"), sep="\t", index=False)

    if len(depths) < min_subjects:
        summary.append(("Status", f"Skipped: fewer than {min_subjects} subjects"))
        finish(pd.DataFrame(columns=RESULT_COLS), pd.DataFrame(columns=["log10_pgen", "log10_c", "phi"]))
        return

    t = d.groupby("junction_aa").agg(k=("subject", "nunique"), pgen=("pgen", "first")).reset_index()
    t["pgen"] = pd.to_numeric(t["pgen"], errors="coerce")
    k_of = t.set_index("junction_aa")["k"]
    anchor_k = anchors["cdr3"].map(k_of).fillna(0).to_numpy()
    anchor_p = pd.to_numeric(anchors["pgen"], errors="coerce").to_numpy()
    ok = np.isfinite(anchor_p) & (anchor_p > 0)
    depth_values = depths.to_numpy(float)
    cal = calibrate(anchor_p[ok], anchor_k[ok], depth_values)
    null = SharingNull(depth_values)

    n_observed = len(t)
    shared = t[(t["k"] >= 2) & (t["pgen"] > 0)].copy()
    p, _ = sharing_pvalues(shared["pgen"].to_numpy(), shared["k"].to_numpy(), null, cal)
    shared["p_value"] = p
    # Singletons take part in the correction: they are hypotheses with p = 1.
    shared["q_value"] = bh_adjust(p, n_observed)
    shared["expected_k"] = expected_sharing(shared["pgen"].to_numpy(), depth_values, cal)
    shared["log10_pgen"] = np.log10(shared["pgen"])
    shared["call"] = np.where(shared["q_value"] < fdr, "selection_candidate", "generation_consistent")
    shared = shared.sort_values(["q_value", "k"], ascending=[True, False])[RESULT_COLS]

    n_called = int((shared["call"] == "selection_candidate").sum())
    summary += [("Public clonotypes (>= 2 subjects, scored)", len(shared)),
                ("Selection candidates", n_called),
                ("Consistent with generation", len(shared) - n_called),
                ("Percent of public clonotypes consistent with generation",
                 round(100 * (len(shared) - n_called) / len(shared), 1) if len(shared) else "NA"),
                ("FDR", fdr), ("Median c", round(float(10 ** np.median(cal.log10_c)), 3)),
                ("Median phi", round(float(np.median(cal.phi)), 3)), ("Status", "Completed")]
    finish(shared, cal.to_frame())


def main():
    ap = argparse.ArgumentParser(description="Test clonotype sharing against a generative null")
    ap.add_argument("--clones", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--subject-col", default="patient")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-subjects", type=int, default=3)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    clones = pd.read_csv(a.clones, sep="\t", usecols=["junction_aa", "sample", "pgen"], low_memory=False)
    run(clones, a.samplesheet, pd.read_csv(a.anchors, sep="\t"), a.subject_col, a.fdr, a.min_subjects, a.outdir)


if __name__ == "__main__":
    main()
