import unittest
import os
from src.gap_filling_dl.biomeneco.model import Model
from bioiso import load
from tests import TEST_DIR


class TestModel(unittest.TestCase):
    def setUp(self):
        # load the test model
        self.model_path = os.path.join(TEST_DIR, "data/models/model_toy_network.xml")
        self.model = load(self.model_path)
        self.model.objective = "e_Biomass__in"

    def test_initialize_model(self):
        model = Model(self.model_path)
        self.assertEqual(len(model.model.reactions), len(self.model.reactions))

    def test_identify_seeds(self):
        # test if identify_seeds method returns a list of tuples
        model = Model(self.model_path)
        self.assertIsInstance(model.seeds, list)
        self.assertIsInstance(model.seeds[0], tuple)

        # test if each tuple in the list contains a string and a string
        for seed in model.seeds:
            self.assertIsInstance(seed[0], str)
            self.assertIsInstance(seed[1], str)

    def test_identify_targets(self):
        model = Model(self.model_path)
        targets = model.identify_targets(objective_function="e_Biomass__in")
        self.assertIsInstance(targets, list)
        if len(targets) > 0:
            self.assertIsInstance(targets[0], tuple)

        for target in targets:
            self.assertIsInstance(target[0], str)
            self.assertIsInstance(target[1], str)

    def test_to_sbml(self):
        try:
            model_path = os.path.join(TEST_DIR, "data/models/model_toy_network.xml")
            model = Model(model_path)
            #print(model.model.reactions)
            model.objective_func = "e_Biomass__in"  # add this line to set the objective function to a string value
            model.to_sbml("test", "tests/data", seeds=False, targets=False)
            # as seeds and targets are False by default, the files should not exist, check if they do not exist
            self.assertFalse(os.path.exists(TEST_DIR + "/data/test_seeds.xml"))  # if the file exists, the test fails
            self.assertFalse(os.path.exists(TEST_DIR + "/data/test_targets.xml"))  # if the file exists, the test fails,
            # however if these lines are commented, the test passes, because it adds an extra "=+1" to file_name

            model.identify_seeds()
            model.identify_targets(objective_function="e_Biomass__in")

            model.to_sbml("test", "/Users/josediogomoura/gap_filling_dl/tests/data", seeds=True, targets=False)
            #print(os.path.exists(TEST_DIR + "/data/test_seeds.xml"))

            model.to_sbml("test", "/Users/josediogomoura/gap_filling_dl/tests/data", seeds=False, targets=True)
            #print(os.path.exists(TEST_DIR + "/data/test_targets.xml"))

            #print(os.path.exists(TEST_DIR + "/data/test_seeds.xml"))
            self.assertTrue(os.path.exists(TEST_DIR + "/data/test_seeds.xml"))

            # if os.path.exists(TEST_DIR + "/data/test_seeds.xml"):
            #     os.remove(TEST_DIR + "/data/test_seeds.xml")
            # if os.path.exists(TEST_DIR + "/data/test_targets.xml"):
            #     os.remove(TEST_DIR + "/data/test_targets.xml")

        except Exception as e:
            print(f"An error occurred: {e}")
            raise e


if __name__ == "__main__":
    unittest.main()
