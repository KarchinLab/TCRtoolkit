"""Manual dev check for tcrdist3_matrix.py's TRBV name conversion.
Not wired into nf-test or CI - run directly: python bin/test_tcrdist3_matrix.py
"""

import unittest
import tcrdist3_matrix as t3m

class TestNameConversion(unittest.TestCase):
    def setUp(self):
        # This method is called before each test
        self.test_cases = [
            ("TCRBV07", "TRBV7*01", "TRBV7*01"),
            ("TCRBV27", "TRBV27*01", "TRBV27*01"),
            ("TCRBV07-02", "TRBV7-2*01", "TRBV7*01"),
            ("TCRBV17-02", "TRBV17-2*01", "TRBV17*01"),
            ("TCRBV10-03*02", "TRBV10-3*02", "TRBV10*02"),
            ("TCRBV07-02*01", "TRBV7-2*01", "TRBV7*01"),
            ("TCRBV10-or09_02*01", "TRBV10/OR9-2*01", "TRBV10/OR9-2*01"),
        ]
        self.split_cases = [
            ("TCRBV06-02/06-03*01", ["TCRBV06-02*01", "TCRBV06-03*01"]),
            ("TCRBV12-03/12-04", ["TCRBV12-03*01","TCRBV12-04*01"])
     ]

    def test_transform_trbv(self):
        for trbv_input, expected_output, _ in self.test_cases:
            result = t3m.transform_trbv(trbv_input)
            self.assertEqual(result, expected_output, f"Failed for input: {trbv_input}")
    
    def test_remove_locus(self):
        for _, trbv_output, remove_locus_output in self.test_cases:
            result = t3m.remove_locus(trbv_output)
            self.assertEqual(result, remove_locus_output, f"Failed for input: {trbv_output}")
    
    def test_split_and_check_genes(self):
        for split_input, split_output in self.split_cases:
            result = t3m.split_and_check_genes(split_input)
            self.assertEqual(result, split_output, f"Failed for input: {split_input}")


class TestClusterMap(unittest.TestCase):
    """cluster_map: clones within the radius are linked, rows sharing a CDR3b are tied, and a
    component needs >= 2 distinct CDR3b to be a cluster."""

    def test_components(self):
        import numpy as np
        import pandas as pd
        import scipy.sparse as sp
        clones = pd.DataFrame({
            "sample":    ["s1", "s1", "s2", "s2", "s1"],
            "cdr3_b_aa": ["CASSA", "CASSB", "CASSA", "CASSZ", "CASSQ"],
            "v_b_gene":  ["TRBV1*01"] * 5,
            "count":     [3, 1, 2, 1, 1],
        })
        # 0-1 within radius; 2 is CASSA again (tied to 0); 3-4 at 30 (beyond 24);
        # a true zero distance is stored as -1 and must count as linked.
        d = np.zeros((5, 5))
        d[0, 1] = d[1, 0] = 12
        d[3, 4] = d[4, 3] = 30
        m = sp.csr_matrix(d)
        out, n_edges = t3m.cluster_map(clones, m, 24, "P1")
        self.assertEqual(n_edges, 1)
        lab = out["tcrdist_cluster"].tolist()
        self.assertEqual(lab[0], lab[1])
        self.assertEqual(lab[0], lab[2])          # same CDR3b, other sample
        self.assertEqual(lab[3], "")              # beyond the radius
        self.assertEqual(lab[4], "")
        self.assertEqual(out.loc[0, "cluster_n_cdr3"], 2)
        self.assertEqual(out.loc[0, "cluster_n_clones"], 3)
        # -1 = zero distance: links 3 and 4 into a cluster.
        d[3, 4] = d[4, 3] = -1
        out, _ = t3m.cluster_map(clones, sp.csr_matrix(d), 24, "P1")
        self.assertNotEqual(out.loc[3, "tcrdist_cluster"], "")
        self.assertEqual(out.loc[3, "tcrdist_cluster"], out.loc[4, "tcrdist_cluster"])


if __name__ == '__main__':
    unittest.main()
