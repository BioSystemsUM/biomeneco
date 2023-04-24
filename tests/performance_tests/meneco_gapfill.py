import os
import unittest
import time
from src.gap_filling_dl.biomeneco.gapfiller import GapFiller
from tests import TEST_DIR


class TestGapFiller(unittest.TestCase):
    def setUp(self):
        self.gf = GapFiller.from_folder(os.path.join(TEST_DIR, "data/6gaps"))
        self.reactions_removed = {
            "model_1": ['R00939_C3_in', 'R03191_C3_in', 'R03066_C3_in', 'R04230_C3_in'],
            "model_2": ['R04440_C3_in', 'R04960_C3_in', 'R01466_C3_in', 'R04954_C3_in', 'R02016_C3_in'],
            "model_3": ['R01231_C3_in', 'R01230_C3_in', 'R01528_C3_in', 'R02749_C3_in'],
            "model_4": ['R01819_C3_in'],
            "model_5": ['R00694_C3_in', 'R00566_C3_in', 'R00228_C3_in', 'R00619_C3_in'],
            "model_6": ['R02739_C3_in', 'R03346_C3_in', 'R05070_C3_in', 'R01863_C3_in', 'R00802_C3_in'],
        }

        self.models_folder = os.path.join(TEST_DIR, "data/original_model")

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

    def test_performance_meneco_manymodels(self):
        for model_name in sorted(os.listdir(self.models_folder)):
            print(f"Running GapFiller for {model_name}")
            model_folder = os.path.join(self.models_folder, model_name)
            gf = GapFiller.from_folder(model_folder)

            start_time = time.time()
            results = gf.run(enumeration=True, json_output=True)
            end_time = time.time()
            elapsed_time = end_time - start_time

            print(f"Time taken for {model_name}: {elapsed_time}")
            print(f"Results for {model_name}:")
            print(results)
            self.print_gapfilling_stats(model_name, results)

    def print_gapfilling_stats(self, model_name, results):
        print(f"--- Gap-filling statistics for {model_name} ---")
        removed_reactions = self.reactions_removed[model_name]
        num_removed_reactions = len(removed_reactions)
        num_reconstructable = len(results["Reconstructable targets"])

        reconstructable_reactions = set(removed_reactions) & set(results["Reconstructable targets"])
        num_reconstructable_removed = len(reconstructable_reactions)

        print(f"Number of removed reactions: {num_removed_reactions}")
        print(f"Number of reconstructable targets: {num_reconstructable}")
        print(f"Number of removed reactions that are reconstructable: {num_reconstructable_removed}")
        print(f"Reconstructable reactions that were removed: {reconstructable_reactions}")
        print("----------------------------------------------------")


if __name__ == '__main__':
    unittest.main()
