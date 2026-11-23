"""Unit tests with known answers for the selection model. Run with: python -m pytest tests"""
import numpy as np
import pytest

import selection_test as sel


def test_poisson_binomial_matches_binomial():
    from scipy import stats
    pmf = sel.poisson_binomial_pmf(np.full((1, 30), 0.2))[0]
    assert np.allclose(pmf, stats.binom.pmf(np.arange(31), 30, 0.2))


def test_poisson_binomial_mixed_probabilities_by_enumeration():
    import itertools
    probs = np.array([0.1, 0.5, 0.9, 0.3])
    want = np.zeros(5)
    for bits in itertools.product([0, 1], repeat=4):
        want[sum(bits)] += np.prod([p if b else 1 - p for p, b in zip(probs, bits)])
    assert np.allclose(sel.poisson_binomial_pmf(probs[None])[0], want)


def test_gamma_nodes_mean_one_and_variance():
    q = sel.gamma_nodes(0.5, 200)
    assert q.mean() == pytest.approx(1)
    assert q.var() == pytest.approx(0.5, rel=0.1)
    assert sel.gamma_nodes(0).tolist() == [1.0]


def _synthetic(c_true=4.0, phi=0.0, n=200000, seed=1):
    rng = np.random.default_rng(seed)
    depths = rng.integers(20000, 80000, 40).astype(float)
    pgen = 10 ** rng.uniform(-9, -5.5, n)
    q = rng.gamma(1 / phi, phi, n) if phi else np.ones(n)
    pi = sel.presence_prob(c_true * q * pgen, depths)
    return pgen, (rng.random(pi.shape) < pi).sum(axis=1), depths, rng


def test_calibration_recovers_constant():
    pgen, k, depths, _ = _synthetic()
    cal = sel.calibrate(pgen, k, depths, fit_phi=False)
    assert np.allclose(10 ** cal.log10_c, 4.0, rtol=0.1)


def test_calibration_recovers_heterogeneity():
    pgen, k, depths, _ = _synthetic(phi=1.0)
    cal = sel.calibrate(pgen, k, depths, trim=1.0)
    assert 0.6 < np.median(cal.phi) < 1.4


@pytest.mark.parametrize("phi", [0.0, 1.0])
def test_pit_is_uniform_on_observed_sequences(phi):
    pgen, k, depths, rng = _synthetic(phi=phi)
    cal = sel.calibrate(pgen[:100000], k[:100000], depths, trim=1.0, fit_phi=phi > 0)
    null = sel.SharingNull(depths)
    seen = k[100000:] >= 1
    u = sel.randomized_pit(pgen[100000:][seen], k[100000:][seen], null, cal, rng)
    assert abs(np.mean(u < 0.05) - 0.05) < 0.01
    assert abs(np.mean(u < 0.5) - 0.5) < 0.02


def test_bh_counts_unscored_hypotheses():
    assert sel.bh_adjust([0.01], m_total=10)[0] == pytest.approx(0.1)
    assert np.allclose(sel.bh_adjust([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])


def test_zero_pgen_is_not_scored():
    depths = np.full(5, 1000.0)
    cal = sel.Calibration(np.array([-10.0, -5.0]), np.zeros(2), np.zeros(2))
    p, _ = sel.sharing_pvalues([0.0, 1e-7], [2, 2], sel.SharingNull(depths), cal)
    assert np.isnan(p[0]) and np.isfinite(p[1])


def test_run_end_to_end(tmp_path):
    import pandas as pd
    rng = np.random.default_rng(3)
    pool = [f"CASS{i:05d}F".replace("0", "A").replace("1", "G") for i in range(4000)]
    pgen = dict(zip(pool, 10 ** rng.uniform(-8, -5, len(pool))))
    rows = [(s, f"s{j}", pgen[s]) for j in range(6) for s in rng.choice(pool, 800, replace=False)]
    rows += [("CASSPUBLICF", f"s{j}", 1e-12) for j in range(6)]
    clones = pd.DataFrame(rows, columns=["junction_aa", "sample", "pgen"])
    sheet = tmp_path / "sheet.csv"
    pd.DataFrame({"sample": [f"s{j}" for j in range(6)], "patient": [f"P{j}" for j in range(6)]}).to_csv(sheet, index=False)
    anchors = pd.DataFrame({"cdr3": pool, "pgen": [pgen[s] for s in pool]})
    sel.run(clones, sheet, anchors, "patient", 0.05, 3, tmp_path)
    out = pd.read_csv(tmp_path / "selection_test.tsv", sep="\t").set_index("junction_aa")
    assert out.loc["CASSPUBLICF", "call"] == "selection_candidate"


def test_too_few_subjects_reports_skip(tmp_path):
    import pandas as pd
    clones = pd.DataFrame({"junction_aa": ["CASSF"] * 2, "sample": ["a", "b"], "pgen": [1e-7] * 2})
    sheet = tmp_path / "s.csv"
    pd.DataFrame({"sample": ["a", "b"], "patient": ["P", "Q"]}).to_csv(sheet, index=False)
    sel.run(clones, sheet, pd.DataFrame({"cdr3": ["CASSF"], "pgen": [1e-7]}), "patient", 0.05, 3, tmp_path)
    s = pd.read_csv(tmp_path / "selection_summary.tsv", sep="\t").set_index("metric")["value"]
    assert s["Status"].startswith("Skipped")


def test_repeat_samples_of_a_patient_count_once(tmp_path):
    import pandas as pd
    clones = pd.DataFrame({"sample": ["a", "b"], "junction_aa": ["X", "X"], "pgen": [1e-7, 1e-7]})
    sheet = tmp_path / "s.csv"
    pd.DataFrame({"sample": ["a", "b"], "patient": ["P", "P"]}).to_csv(sheet, index=False)
    d, unit = sel.subject_table(clones, sheet, "patient")
    assert unit == "patient" and d["subject"].nunique() == 1
