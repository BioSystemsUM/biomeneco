import unittest
import os
from src.gap_filling_dl.biomeneco.model import Model
from bioiso import load


class TestModel(unittest.TestCase):
    def setUp(self):
        # load the test model
        self.model_path = "/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/models/toy_network.xml"
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
        model = Model(self.model_path)
        model.to_sbml("test.xml", "tests/performance_tests/data/models")
        self.assertTrue(os.path.exists("tests/performance_tests/data/models/test.xml"))
        os.remove("tests/performance_tests/data/models/test.xml")


if __name__ == "__main__":
    unittest.main()
