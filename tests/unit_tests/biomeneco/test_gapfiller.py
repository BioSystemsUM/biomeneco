import unittest
import os
from unittest.mock import patch, MagicMock, mock_open
import io
import cobra

from src.gap_filling_dl.biomeneco.gapfiller import GapFiller
import unittest.mock as mock


class TestGapFiller(unittest.TestCase):

    # def __init__(self, methodName: str = ...):
    #     super().__init__(methodName)
    #     self.test_folder =

    def setUp(self):
        # Set up a test folder containing example files for gap-filling
        self.model_path = "/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/6gaps_model/model.xml"
        self.universal_model_path = "/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/6gaps_model/universal_model.xml"
        self.gf = GapFiller(self.model_path, self.universal_model_path)

        self.seeds_path = "/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/6gaps_model/seeds.xml"
        self.targets_path = "/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/6gaps_model/targets.xml"

        self.folder_path = "/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/6gaps_model"

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
            gf = GapFiller(self.model_path, self.universal_model_path)
            self.assertEqual(gf.seeds_path, '/path/to/seeds')

    def test_targets_path(self):
        with mock.patch.object(GapFiller, 'targets_path', '/path/to/targets'):
            gf = GapFiller(self.model_path, self.universal_model_path)
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
        self.gf.run_meneco = MagicMock()

        expected_result = {'Draft network file': mock.ANY,
                           'Seeds file': mock.ANY,
                           'Targets file': mock.ANY,
                           'Unproducible targets': mock.ANY}

        self.gf.run_meneco.return_value = expected_result

        print(self.gf.run_meneco.__dict__)

        result = self.gf.run(enumeration=True, json_output=True)
        self.gf.run_meneco.assert_called_once_with(enumeration=True, json_output=True)
        self.assertEqual(result, expected_result)

    def test_run_meneco(self):
        self.gf._seeds_path = self.seeds_path
        self.gf._targets_path = self.targets_path
        print(self.gf._seeds_path)
        print(self.gf._targets_path)
        expected_result = {'Draft network file': mock.ANY,
                           'Seeds file': mock.ANY,
                           'Targets file': mock.ANY,
                           'Unproducible targets': mock.ANY}

        result = self.gf.run_meneco(enumeration=True, json_output=True)
        self.assertEqual(result, expected_result)

    def test_from_folder(self):
        # folder path is: /Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/6gaps_model

        model_file = "model.xml"
        seeds_file = "seeds.xml"
        targets_file = "targets.xml"
        universal_model_file = "universal_model.xml"
        for file in [model_file, seeds_file, targets_file, universal_model_file]:
            open(os.path.join(self.folder_path, file), 'w').close()

        gap_filler = GapFiller.from_folder(self.folder_path)

        self.assertIsNotNone(gap_filler)
        self.assertEqual(gap_filler.seeds_path, os.path.join(self.folder_path, seeds_file))
        self.assertEqual(gap_filler.targets_path, os.path.join(self.folder_path, targets_file))

    def test_from_folder_missing_files(self):
        # Ensure that the required files are not present
        model_file = "model.xml"
        seeds_file = "seeds.xml"
        targets_file = "targets.xml"
        universal_model_file = "universal_model.xml"
        missing_files = [model_file, seeds_file, targets_file, universal_model_file]
        existing_files = os.listdir(self.folder_path)
        for file in existing_files:
            if file in missing_files:
                os.remove(os.path.join(self.folder_path, file))

        # Ensure that the method raises FileNotFoundError
        with self.assertRaises(FileNotFoundError):
            GapFiller.from_folder(self.folder_path)

    def test_evalute_results(self):
        pass


if __name__ == '__main__':
    unittest.main()
