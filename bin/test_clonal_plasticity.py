import numpy as np
import pandas as pd
import pytest

import clonal_plasticity as pl


def test_jsd_bounds():
    assert pl.jsd([1, 0], [1, 0]) == pytest.approx(0)
    assert pl.jsd([1, 0], [0, 1]) == pytest.approx(1)
    assert pl.jsd([3, 1], [1, 3]) == pytest.approx(pl.jsd([1, 3], [3, 1]))


def _cells(n_clones, n, p1, p2, rng, background=None):
    rows = []
    for c in range(n_clones):
        for t, p in (("d0", p1), ("d30", p2)):
            for s in rng.choice(["A", "B"], n, p=[p, 1 - p]):
                rows.append((f"c{c}", t, s))
    return pd.DataFrame(rows, columns=["clone", "timepoint", "state"])


def test_stable_clones_give_uniform_time_pvalues():
    rng = np.random.default_rng(1)
    out = pl.clone_plasticity(_cells(300, 15, 0.5, 0.5, rng), n_boot=50, n_perm=199, seed=2)
    # JSD on two states and 15 cells takes few values, so the test is valid but conservative.
    assert (out["p_time"] < 0.05).mean() <= 0.08


def test_shifted_clone_is_detected_with_interval_covering_truth():
    rng = np.random.default_rng(1)
    cells = _cells(1, 200, 0.9, 0.1, rng)
    out = pl.clone_plasticity(cells, n_boot=500, n_perm=199, seed=3).iloc[0]
    truth = float(pl.jsd([0.9, 0.1], [0.1, 0.9]))
    assert out["p_time"] < 0.01 and out["ci_low"] <= truth <= out["ci_high"]


def test_state_null_discounts_population_wide_shift():
    rng = np.random.default_rng(4)
    # Every clone follows the same population shift: real change, but nothing clone-specific.
    out = pl.clone_plasticity(_cells(60, 20, 0.8, 0.2, rng), n_boot=50, n_perm=199, seed=5)
    assert (out["p_time"] < 0.05).mean() > 0.9
    assert (out["p_state"] < 0.05).mean() < 0.15


def test_clone_below_min_cells_is_skipped_and_order_is_respected():
    cells = pd.DataFrame({"clone": ["a"] * 8, "timepoint": ["d0"] * 6 + ["d9"] * 2, "state": ["A"] * 8})
    assert pl.clone_plasticity(cells, min_cells=5).empty
    two = pd.DataFrame({"clone": ["a"] * 12, "timepoint": ["late"] * 6 + ["early"] * 6,
                        "state": ["A", "B"] * 6})
    row = pl.clone_plasticity(two, order=["early", "late"], n_boot=20, n_perm=20).iloc[0]
    assert (row["t1"], row["t2"]) == ("early", "late")


def test_run_takes_timepoints_from_samplesheet(tmp_path):
    cells = pd.DataFrame({"sample": ["s1"] * 6 + ["s2"] * 6, "patient": "P", "clone_id": "c", "annot": ["A"] * 6 + ["B"] * 6})
    sheet = pd.DataFrame({"sample": ["s1", "s2"], "timepoint": ["d0", "d30"], "origin": ["tumor", "tumor"]})
    pl.run(cells, sheet, "d0,d30", 5, 50, 50, 1, tmp_path)
    out = pd.read_csv(tmp_path / "clonal_plasticity.tsv", sep="\t").iloc[0]
    assert out["t1"] == "d0" and out["jsd"] == pytest.approx(1)


def test_no_eligible_clone_writes_header(tmp_path):
    cells = pd.DataFrame({"sample": ["s1"], "patient": "P", "clone_id": "c", "annot": "A"})
    sheet = pd.DataFrame({"sample": ["s1"], "timepoint": ["d0"], "origin": ["tumor"]})
    pl.run(cells, sheet, "", 5, 10, 10, 1, tmp_path)
    assert "jsd" in (tmp_path / "clonal_plasticity.tsv").read_text()


def test_numeric_timepoints_sort_numerically():
    assert pl.timepoint_order(["10", "2"], "") == ["2", "10"]
    assert pl.timepoint_order(["a", "b"], "b,a") == ["b", "a"]
