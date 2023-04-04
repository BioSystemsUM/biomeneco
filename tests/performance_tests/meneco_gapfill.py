import os
import unittest
import time
from src.gap_filling_dl.biomeneco.gapfiller import GapFiller
from tests import TEST_DIR


class TestGapFiller(unittest.TestCase):
    def setUp(self):
        self.gf = GapFiller.from_folder(os.path.join(TEST_DIR, "data/6gaps"))

    def test_performance_meneco(self):
        start_time = time.time()
        results = self.gf.run(enumeration=True, json_output=True)
        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f"time taken: {elapsed_time}")
        print(results)
        self.assertTrue(isinstance(results, dict))
        self.assertTrue("Draft network file" in results)
        self.assertTrue("Seeds file" in results)
        self.assertTrue("Targets file" in results)
        self.assertTrue("Unproducible targets" in results)
        self.assertTrue("Repair db file" in results)
        self.assertTrue("Unreconstructable targets" in results)
        self.assertTrue("Reconstructable targets" in results)
        self.assertTrue(isinstance(results["Unproducible targets"], list))
        self.assertTrue(isinstance(results["Unreconstructable targets"], list))
        self.assertTrue(isinstance(results["Reconstructable targets"], list))


