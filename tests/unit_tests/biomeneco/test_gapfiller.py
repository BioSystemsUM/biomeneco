import unittest
import os
from unittest.mock import patch, MagicMock, mock_open
import io
import cobra
from src.gap_filling_dl.biomeneco.gapfiller import GapFiller
import unittest.mock as mock

from tests import TEST_DIR


class TestGapFiller(unittest.TestCase):

    def setUp(self):
        # Set up a test folder containing example files for gap-filling
        self.model_path = os.path.join(TEST_DIR, "data/test_model/model_toy_network.xml")

        self.universal_model_path = os.path.join(TEST_DIR, "data/test_model/universal_model_kegg.xml")
        self.gf = GapFiller(self.model_path, self.universal_model_path)

        self.seeds_path = os.path.join(TEST_DIR, "data/test_model/test_seeds.xml")
        self.targets_path = os.path.join(TEST_DIR, "data/test_model/test_targets.xml")

        self.folder_path = os.path.join(TEST_DIR, "data/test_model")

    # def tearDown(self):
    #     # remove the mock files
    #     os.remove(self.model_path)
    #     os.remove(self.universal_model_path)

    def test_init(self):
        self.assertEqual(self.gf.model_path, self.model_path)
        self.assertEqual(self.gf.universal_model_path, self.universal_model_path)
        self.assertIsNone(self.gf.results_meneco)
        self.assertIsNone(self.gf._seeds_path)
        self.assertIsNone(self.gf._targets_path)

    def test_seeds_path(self):
        with mock.patch.object(GapFiller, 'seeds_path', '/path/to/seeds'):
            gf = GapFiller('/path/to/model', '/path/to/universal_model')
            self.assertEqual(gf.seeds_path, '/path/to/seeds')

    def test_targets_path(self):
        with mock.patch.object(GapFiller, 'targets_path', '/path/to/targets'):
            gf = GapFiller('/path/to/model', '/path/to/universal_model')
            self.assertEqual(gf.targets_path, '/path/to/targets')

    @mock.patch('cobra.io.read_sbml_model')  # Mock the read_sbml_model function from cobra
    def test_universal_model(self, mock_read_sbml_model):
        mock_model = mock.Mock(spec=cobra.Model)
        mock_read_sbml_model.return_value = mock_model

        result = self.gf.universal_model
        mock_read_sbml_model.assert_called_once_with(
            self.gf.universal_model_path)  # Check that the correct file is read

        self.assertEqual(result, mock_model)  # Check that the correct model is returned

    def test_run(self):
        with mock.patch.object(GapFiller, 'run') as mock_run_meneco:
            self.gf.model_path = self.model_path
            self.gf.universal_model_path = self.universal_model_path
            self.gf._seeds_path = self.seeds_path
            self.gf._targets_path = self.targets_path

            print(self.gf.model_path)
            print(self.gf.universal_model_path)
            print(self.gf._seeds_path)
            print(self.gf._targets_path)

            result = self.gf.run(enumeration=True, json_output=True)
            mock_run_meneco.assert_called_once_with(enumeration=True, json_output=True)

            self.assertEqual(result, mock_run_meneco.return_value)

    def test_run_meneco(self):
        # similar as test_run but with the run_meneco method...
        with mock.patch.object(GapFiller, 'run_meneco') as mock_run_meneco:
            self.gf.model_path = self.model_path
            self.gf.universal_model_path = self.universal_model_path
            self.gf._seeds_path = self.seeds_path
            self.gf._targets_path = self.targets_path

            result = self.gf.run_meneco(enumeration=True, json_output=True)
            mock_run_meneco.assert_called_once_with(enumeration=True, json_output=True)

            self.assertEqual(result, mock_run_meneco.return_value)

    def test_from_folder(self):
        folder_path = os.path.join(TEST_DIR, "data/test_model")

        model_file = "model_toy_network.xml"
        seeds_file = "test_seeds.xml"
        targets_file = "test_targets.xml"
        universal_model_file = "universal_model_kegg.xml"

        # check if the files are present
        # Check if all required files exist in the folder
        for file in [model_file, seeds_file, targets_file, universal_model_file]:
            file_path = os.path.join(folder_path, file)
            self.assertTrue(os.path.isfile(file_path),
                            f"File {file_path} does not exist in folder {folder_path}")

        gap_filler = GapFiller.from_folder(folder_path)

        # Check if the correct files are assigned to the gap-filler
        self.assertIsNotNone(gap_filler.model_path)
        self.assertIsNotNone(gap_filler.universal_model_path)
        self.assertIsNotNone(gap_filler._seeds_path)
        self.assertIsNotNone(gap_filler._targets_path)

        # Check if the correct files are assigned to the gap-filler
        self.assertEqual(gap_filler.model_path, os.path.join(folder_path, model_file))
        self.assertEqual(gap_filler.universal_model_path, os.path.join(folder_path, universal_model_file))
        self.assertEqual(gap_filler._seeds_path, os.path.join(folder_path, seeds_file))
        self.assertEqual(gap_filler._targets_path, os.path.join(folder_path, targets_file))

    def test_evalute_results(self):
        gap_filler = GapFiller(self.model_path, self.universal_model_path)
        gap_filler._seeds_path = self.seeds_path
        gap_filler._targets_path = self.targets_path

        results_meneco = gap_filler.run(enumeration=True, json_output=True)
        print(results_meneco)

        # Check output is a dictionary
        self.assertIsInstance(results_meneco, dict)

        # Check that the output dictionary has the expected keys
        expected_keys = {'Draft network file', 'Seeds file', 'Targets file', 'Unproducible targets'}
        self.assertSetEqual(set(results_meneco.keys()), expected_keys)

        evaluted_results = gap_filler.evaluate_results()
        expected_evaluated_keys = {'Seed reactions', 'Target reactions', 'Draft reactions', 'Gap-filled reactions',
                                   'Unproducible reactions'}

        # Check that the output dictionary has the expected keys
        self.assertSetEqual(set(evaluted_results.keys()), expected_evaluated_keys)

        # test if verbose is True, print the results
        with mock.patch('builtins.print') as mock_print:
            gap_filler.evaluate_results(verbose=True)
            mock_print.assert_called()


if __name__ == '__main__':
    unittest.main()
