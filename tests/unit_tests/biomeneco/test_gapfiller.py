import unittest
import os
from unittest.mock import patch, MagicMock, mock_open
import io
import cobra
from src.gap_filling_dl.biomeneco.gapfiller import GapFiller
import unittest.mock as mock
from cobra import Reaction

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

        self.model_copy_path = os.path.join(TEST_DIR, "data/test_model/model_toy_network_copy.xml")

        # remove reactions frmo the copy
        self.model_copy = cobra.io.read_sbml_model(self.model_copy_path)

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

    # def test_add_reactions_to_model(self):
    #     # reactions to remove:
    #     reactions_to_remove = ["R00200__in"]  # EX_C00267__in does not exist in the universal model
    #     copy_model = cobra.io.read_sbml_model(self.model_copy_path)
    #     copy_model.objective = 'e_Biomass__in'
    #     # remove permenantly the reactions from the model:
    #     for r_id in reactions_to_remove:
    #         copy_model.reactions.get_by_id(r_id).remove_from_model()
    #
    #     # check if the reactions are in our model 'EX_C00267__in', 'R01061__in', 'R01015__in'
    #     for r_id in reactions_to_remove:
    #         self.assertNotIn(r_id, copy_model.reactions)
    #
    #     # add the reactions to the model
    #     gap_filler = GapFiller(self.model_copy_path, self.universal_model_path)
    #
    #     gap_filler.add_reactions_to_model(reactions_to_remove)
    #
    #     # check if the reactions are in our model 'EX_C00267__in', 'R01061__in', 'R01015__in'
    #     for r_id in reactions_to_remove:
    #         self.assertIn(r_id, copy_model.reactions)

    def test_add_reactions_to_model(self):
        # reactions to add:
        reactions_to_add = ["R00002__in", "R00010__in", "R00019__in"]
        model = cobra.io.read_sbml_model(self.model_copy_path)
        model.objective = 'e_Biomass__in'

        # add the reactions to the model
        gap_filler = GapFiller(self.model_copy_path, self.universal_model_path)
        gap_filler.add_reactions_to_model(reactions_to_add)

        # check that the reactions were added to the model
        for reaction_id in reactions_to_add:
            self.assertIn(reaction_id, [r.id for r in model.reactions])

        # check that the cobra_filled_model object is not empty
        self.assertIsNotNone(gap_filler.cobra_filled_model)

    def test_compare_reactions(self):
        # create initial model reactions
        initial_reactions = ["R1", "R2", "R3"]

        # create filled model reactions
        filled_model_reactions = ["R1", "R2", "R3", "R4", "R5"]

        # create a mock for the filled model
        cobra_filled_model_mock = mock.MagicMock()
        cobra_filled_model_mock.reactions = [mock.MagicMock(id=id) for id in filled_model_reactions]

        # create a mock for the cobra.io.read_sbml_model function
        read_sbml_model_mock = mock.MagicMock(
            return_value=mock.MagicMock(reactions=[mock.MagicMock(id=id) for id in initial_reactions]))

        # create GapFiller object with mocks
        with mock.patch('cobra.io.read_sbml_model', read_sbml_model_mock):
            gap_filler = GapFiller(model_path=self.model_path, universal_model_path=self.universal_model_path)
            gap_filler.cobra_filled_model = cobra_filled_model_mock

            # test added reactions
            added_reactions = gap_filler.compare_reactions()
            expected_added_reactions = ["R4", "R5"]

            # assert that added reactions match expected added reactions
            self.assertCountEqual(added_reactions, expected_added_reactions)

    def test_compare_reactions_no_added_reactions(self):
        # Create a mock filled model with no additional reactions
        filled_model = cobra.io.read_sbml_model(self.model_path)
        self.gf.cobra_filled_model = filled_model

        # Call compare_reactions and check the returned value
        added_reactions = self.gf.compare_reactions()
        self.assertListEqual(added_reactions, [])


if __name__ == '__main__':
    unittest.main()
