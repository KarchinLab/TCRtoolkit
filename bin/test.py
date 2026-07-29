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


if __name__ == '__main__':
    unittest.main()
