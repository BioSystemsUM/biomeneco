import glob
import tempfile
import unittest
import os

import cobra.io

from src.gap_filling_dl.biomeneco.model import Model
from bioiso import load
from tests import TEST_DIR


class TestModel(unittest.TestCase):
    def setUp(self):
        # load the test model
        self.model_path = os.path.join(TEST_DIR, "data/original_model/iDS372.xml")
        self.model = Model.from_sbml(self.model_path, objective_function="Biomass_assembly_C3_in")
        self.parent_folder = '/Users/josediogomoura/gap_filling_dl/tests/data/original_model'
        self.universal_model_path = "/Users/josediogomoura/gap_filling_dl/tests/data/original_model/universal_model_kegg.xml"

    def test_from_sbml(self):
        model = Model.from_sbml(self.model_path, objective_function="Biomass_assembly_C3_in")
        self.assertEqual(len(model.reactions), len(self.model.reactions))

    def test_identify_seeds(self):
        # check that the output is a list
        self.assertIsInstance(self.model.identify_seeds(), list)
        # check that the list contains tuples of length 2
        for seed in self.model.identify_seeds():
            self.assertIsInstance(seed, tuple)
            self.assertEqual(len(seed), 2)
        # check that the list contains the expected number of seeds
        self.assertEqual(len(self.model.identify_seeds()), 67)

    def test_identify_targets(self):
        # check that the output is a list
        self.assertIsInstance(self.model.identify_targets(), list)
        # check that the list contains tuples of length 2
        for target in self.model.identify_targets():
            self.assertIsInstance(target, tuple)
            self.assertEqual(len(target), 2)
        # check that the list contains the expected number of targets
        self.assertEqual(len(self.model.identify_targets()), 0)

    def test_identify_targets_removed_reactions(self):
        # make a copy of the model
        sixgaps_model = '/Users/josediogomoura/gap_filling_dl/tests/data/6gaps/6gaps_model.xml'
        model = Model.from_sbml(sixgaps_model, objective_function="Biomass_assembly_C3_in")

        # run identify targets
        targets = model.identify_targets()
        #print(targets)

        # check that the output is a list
        self.assertIsInstance(targets, list)
        # check that the list contains tuples of length 2
        for target in targets:
            self.assertIsInstance(target, tuple)
            self.assertEqual(len(target), 2)
        # check that the list contains the expected number of targets
        self.assertEqual(len(targets), 7)

    def test_to_sbml(self):
        # create temporary folder for testing
        with tempfile.TemporaryDirectory() as tmpdir:
            # test saving seeds and targets files
            self.model.to_sbml("test_model.xml", save_path=tmpdir, seeds=True, targets=True)
            self.assertTrue(os.path.exists(os.path.join(tmpdir, "test_model_seeds.xml")))
            self.assertTrue(os.path.exists(os.path.join(tmpdir, "test_model_targets.xml")))
            # test not saving seeds and targets files
            self.model.to_sbml("test_model_1.xml", save_path=tmpdir, seeds=False, targets=False)
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "test_model_1_seeds.xml")))
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "test_model_1_targets.xml")))

    def test_create_random_knockout_models_generates_models(self):
        num_models = 6
        min_knockouts = 1
        max_knockouts = 6

        # print(f"Reactions: {self.model.reactions}")2

        knockout_models = self.model.create_random_knockout_models(num_models, min_knockouts, max_knockouts,
                                                                   universal_model=self.universal_model_path)
        print(knockout_models)
        self.assertEqual(len(knockout_models), num_models)

    def test_create_random_knockout_models_generates_seeds_targets(self):
        num_models = 6
        min_knockouts = 1
        max_knockouts = 6
        self.model.create_random_knockout_models(num_models, min_knockouts, max_knockouts, generate_files=True,
                                                 parent_folder=self.parent_folder, seeds_targets=True,
                                                 universal_model=self.universal_model_path)
        for i in range(1, num_models + 1):
            model_folder = os.path.join(self.parent_folder, f"model_{i}")
            seeds_files_pattern = os.path.join(model_folder, f"{self.model.id}_knockout_*_seeds.xml")
            targets_files_pattern = os.path.join(model_folder, f"{self.model.id}_knockout_*_targets.xml")

            seeds_files = glob.glob(seeds_files_pattern)
            targets_files = glob.glob(targets_files_pattern)

            self.assertTrue(len(seeds_files) > 0, f"No seed files found in folder {model_folder}")
            self.assertTrue(len(targets_files) > 0, f"No target files found in folder {model_folder}")


if __name__ == "__main__":
    unittest.main()
