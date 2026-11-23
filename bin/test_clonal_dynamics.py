"""A clone with several nucleotide rows in a sample must be summed, not averaged."""
import numpy as np
import pandas as pd

import clonal_dynamics as cd


def _clonotypes():
    rows = []
    for timepoint, rank, clones in (("pre", 1, [("A", 3), ("A", 2), ("B", 50)]),
                                    ("post", 2, [("A", 40), ("B", 50)])):
        for cdr3, n in clones:
            rows.append({"junction_aa": cdr3, "v_call": "V", "j_call": "J", "duplicate_count": n,
                         "patient": "P", "origin": "tumor", "timepoint": timepoint, "timepoint_rank": rank})
    d = pd.DataFrame(rows)
    d["frequency"] = d["duplicate_count"] / d.groupby("timepoint")["duplicate_count"].transform("sum")
    return d


def test_duplicate_rows_are_summed_and_q_values_survive():
    _, detail, _ = cd.compute_dynamics(_clonotypes(), "patient", "origin", "timepoint", "up", 2.0, 0.05)
    a = detail[detail["junction_aa"] == "A"].iloc[0]
    assert a["count_pre"] == 5 and a["count_pre"] % 1 == 0
    assert np.isfinite(a["p_value"]) and np.isfinite(a["q_value"])
