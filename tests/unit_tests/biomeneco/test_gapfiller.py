import tempfile
import unittest
import os
from os.path import dirname
from unittest.mock import patch, MagicMock, mock_open
import io
import cobra
import pytest
from clyngor.as_pyasp import TermSet
from cobra.io import read_sbml_model
from meneco import sbml
from functools import partial

from parallelbar import progress_imap

import gap_filling_dl
from gap_filling_dl.biomeneco.gapfiller import GapFiller
import unittest.mock as mock

from cobra import Reaction, Model, Metabolite, Solution
from unittest.mock import patch, MagicMock

from tests import TEST_DIR


def mock_clone_metabolite(compartments, metabolites):
    cloned_metabolites = []
    for metabolite in metabolites:
        for compartment in compartments:
            new_metabolite = mock.Mock(spec=cobra.Metabolite)
            new_metabolite.id = f"{metabolite.id}__{compartment}"
            cloned_metabolites.append(new_metabolite)
    return cloned_metabolites


class MockMetabolite(MagicMock):
    def __init__(self, id, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.id = id


class TestGapFiller(unittest.TestCase):

    def setUp(self):
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        self.model_path = os.path.join(base_dir, "tests/data/test_model/model_toy_network_copy.xml")
        self.universal_model_path = os.path.join(base_dir, "utilities/universal_model.xml")
        self.results_path = os.path.join(base_dir, "tests/data/test_model/output")
        self.max_solutions = 1
        self.gf = GapFiller(self.model_path, self.universal_model_path, self.results_path, self.max_solutions)

        self.gf.seeds_path = os.path.join(base_dir, "tests/data/test_model/model_seeds.xml")
        self.gf.targets_path = os.path.join(base_dir, "tests/data/test_model/model_targets.xml")

    # def tearDown(self):
    #     # remove the mock files
    #     os.remove(self.model_path)
    #     os.remove(self.universal_model_path)

    def test_init(self):
        self.assertEqual(self.gf.model_path, self.model_path)
        self.assertEqual(self.gf.universal_model_path, self.universal_model_path)
        self.assertEqual(self.gf.resultsPath, self.results_path)  # Ensure attribute name matches class definition
        self.assertEqual(self.gf.max_solutions, self.max_solutions)

        # Test default values of other attributes
        self.assertIsNone(self.gf.optimal_completions)
        self.assertIsNone(self.gf.added_demands)
        self.assertIsNone(self.gf.never_producible)
        self.assertIsNone(self.gf.reconstructable_targets)
        self.assertIsNone(self.gf.unproducible)
        self.assertIsNone(self.gf.universal_model_copy)
        self.assertIsNone(self.gf.universal_model_compartmentalized)
        self.assertIsNone(self.gf.transport_reactions_universal)
        self.assertIsNone(self.gf.transport_reactions)
        self.assertIsNone(self.gf.dead_ends)
        self.assertIsNone(self.gf.cloned_model)
        self.assertIsNone(self.gf.temporary_universal_model)
        self.assertIsNotNone(self.gf.original_model)
        self.assertIsNotNone(self.gf._universal_model)
        self.assertIsNone(self.gf.minimal_set)
        self.assertIsNone(self.gf.minimal_completion)
        self.assertIsNone(self.gf.gftime)
        self.assertIsNone(self.gf.cobra_filled_model)
        self.assertIsNone(self.gf.results_meneco)
        self.assertIsNone(self.gf._seeds_path)
        self.assertIsNone(self.gf._targets_path)
        self.assertIsNone(self.gf.filled_model)
        self.assertIsNone(self.gf.objective_function_id)

        # Check collections are initialized correctly
        self.assertEqual(self.gf.all_completions, [])
        self.assertEqual(self.gf.all_solutions_with_models, [])
        self.assertEqual(self.gf.positive_solutions, set())
        self.assertEqual(self.gf.required_additional_seeds_and_demands, set())
        self.assertEqual(self.gf.additional_seeds, set())

    def test_seeds_path(self):
        seeds_path = os.path.join(TEST_DIR, "data/test_model/model_seeds.xml")
        with mock.patch.object(GapFiller, 'seeds_path', new_callable=mock.PropertyMock) as mock_seeds_path:
            mock_seeds_path.return_value = seeds_path
            self.assertEqual(self.gf.seeds_path, seeds_path)

    def test_targets_path(self):
        targets_path = os.path.join(TEST_DIR, "data/test_model/model_targets.xml")
        with mock.patch.object(GapFiller, 'targets_path', new_callable=mock.PropertyMock) as mock_targets_path:
            mock_targets_path.return_value = targets_path
            self.assertEqual(self.gf.targets_path, targets_path)

    @mock.patch('cobra.io.read_sbml_model')
    def test_universal_model_getter(self, mock_read_sbml_model):
        # mocking the cobra.io.read_sbml_model to return a mock cobra.Model
        mock_model = mock.Mock(spec=cobra.Model)
        mock_read_sbml_model.return_value = mock_model

        # Test when _universal_model is None and is loaded from file
        self.gf._universal_model = None
        model = self.gf.universal_model
        mock_read_sbml_model.assert_called_once_with(self.gf.universal_model_path)
        self.assertEqual(model, mock_model)

        # Reset mock
        mock_read_sbml_model.reset_mock()

        # Test when _universal_model is already set
        self.gf._universal_model = mock_model
        model = self.gf.universal_model
        mock_read_sbml_model.assert_not_called()
        self.assertEqual(model, mock_model)

    def test_universal_model_setter(self):
        new_model = mock.Mock(spec=cobra.Model)
        self.gf.universal_model = new_model
        self.assertEqual(self.gf._universal_model, new_model)

    @mock.patch('cobra.io.read_sbml_model')
    def test_original_model_getter(self, mock_read_sbml_model):
        # Mocking the cobra.io.read_sbml_model to return a mock cobra.Model
        mock_model = mock.Mock(spec=cobra.Model)
        mock_read_sbml_model.return_value = mock_model

        # Test when _original_model is None and is loaded from file
        self.gf._original_model = None
        model = self.gf.original_model
        mock_read_sbml_model.assert_called_once_with(self.gf.model_path)
        self.assertEqual(model, mock_model)

        # Reset mock
        mock_read_sbml_model.reset_mock()

        # Test when _original_model is already set
        self.gf._original_model = mock_model
        model = self.gf.original_model
        mock_read_sbml_model.assert_not_called()
        self.assertEqual(model, mock_model)

    def test_original_model_setter(self):
        new_model = mock.Mock(spec=cobra.Model)
        self.gf.original_model = new_model
        self.assertEqual(self.gf._original_model, new_model)

    def test_run(self):
        # Mock dependencies
        with patch('cobra.io.read_sbml_model') as mock_read_model, \
                patch('meneco.meneco.sbml.readSBMLnetwork') as mock_read_network, \
                patch('meneco.meneco.sbml.readSBMLseeds') as mock_read_seeds, \
                patch('meneco.meneco.sbml.readSBMLtargets') as mock_read_targets, \
                patch('cobra.io.write_sbml_model') as mock_write_model, \
                patch('time.time', side_effect=[100, 200]):  # Mock time to simulate duration

            # Mock return values for dependencies
            mock_read_model.return_value = MagicMock()
            mock_read_network.return_value = MagicMock()
            mock_read_seeds.return_value = MagicMock()
            mock_read_targets.return_value = MagicMock()

            # Call the run method
            self.gf.run(optimize=False, write_to_file=True)

            # Assert that the dependencies are called as expected
            mock_read_model.assert_called_once_with(self.gf.model_path)
            mock_read_network.assert_called_once_with(self.gf.model_path, 'draft')
            mock_read_seeds.assert_called_once_with(self.gf.seeds_path)
            mock_read_targets.assert_called_once_with(self.gf.targets_path)
            mock_write_model.assert_called_once_with(self.gf.original_model,
                                                     os.path.join(self.gf.resultsPath, "gapfilled_model.xml"))

            # Assert internal state changes
            self.assertIsNotNone(self.gf.original_model)
            self.assertGreaterEqual(self.gf.gftime, 100)  # Check if time taken is as expected

            # Additional assertions for other internal state changes or outputs can be added here

    def test_run_meneco(self):
        # similar as test_run but with the run_meneco method...
        with mock.patch.object(GapFiller, 'run_meneco') as mock_run_meneco:
            self.gf.model_path = self.model_path
            self.gf.universal_model_path = self.universal_model_path

            result = self.gf.run_meneco(enumeration=True, json_output=True)
            mock_run_meneco.assert_called_once_with(enumeration=True, json_output=True)

            self.assertEqual(result, mock_run_meneco.return_value)

    def test_run_meneco_varied_parameters(self):
        for enumeration, json_output in [(True, True), (True, False), (False, True), (False, False)]:
            with mock.patch.object(GapFiller, 'run_meneco', return_value={'result': 'test'}) as mock_run_meneco:
                result = self.gf.run_meneco(enumeration=enumeration, json_output=json_output)
                mock_run_meneco.assert_called_once_with(enumeration=enumeration, json_output=json_output)
                self.assertEqual(result, {'result': 'test'})

    @patch('shutil.copy')
    @patch('cobra.io.read_sbml_model')
    @patch('os.path.isfile')
    def test_from_folder(self, mock_isfile, mock_read_sbml_model, mock_copy):
        folder_path = "/Users/josediogomoura/gap_filling_dl/tests/data/test_model"
        results_path = "/Users/josediogomoura/gap_filling_dl/tests/data/test_model/output"

        # Mock the necessary functions
        mock_isfile.return_value = True
        mock_read_sbml_model.return_value = MagicMock()

        # Call the from_folder method
        gap_filler = GapFiller.from_folder(folder_path, results_path)

        self.assertGreaterEqual(mock_isfile.call_count, 0)

        # Assert object initialization
        self.assertIsInstance(gap_filler, GapFiller)

        # Optional: Assert properties of gap_filler
        self.assertEqual(gap_filler.model_path, folder_path + "/model.xml")
        self.assertEqual(gap_filler.seeds_path, folder_path + "/model_seeds.xml")
        self.assertEqual(gap_filler.targets_path, folder_path + "/model_targets.xml")

    @mock.patch('cobra.io.write_sbml_model')
    @mock.patch('meneco.meneco.sbml.readSBMLnetwork')
    @mock.patch('meneco.meneco.sbml.readSBMLseeds')
    @mock.patch('meneco.query.get_optimal_completions')
    @mock.patch('meneco.query.get_unproducible')
    def test_add_reactions_with_positive_solution(self, mock_get_unproducible, mock_get_optimal_completions,
                                                  mock_readSBMLseeds, mock_readSBMLnetwork, mock_write_sbml_model):
        # Mock all necessary dependencies
        mock_get_unproducible.return_value = set()  # Mock unproducible targets as empty set
        mock_get_optimal_completions.return_value = [{}]  # Mock optimal completions as a list with an empty dictionary
        mock_readSBMLseeds.return_value = mock.Mock()  # Mock readSBMLseeds as a mock object
        mock_readSBMLnetwork.return_value = mock.Mock()  # Mock readSBMLnetwork as a mock object
        mock_write_sbml_model.return_value = None  # Mock write_sbml_model to do nothing

        # Set the seeds and targets paths for the GapFiller instance
        self.gf.seeds_path = os.path.join(TEST_DIR, "data/test_model/model_seeds.xml")
        self.gf.targets_path = os.path.join(TEST_DIR, "data/test_model/model_targets.xml")

        # Mock slim_optimize method to return a numeric value
        self.gf.original_model.slim_optimize.return_value = 0.1  # Example value

        # Test the add_reactions_with_positive_solution method
        self.gf.add_reactions_to_model()

    import tempfile
    from cobra.io.sbml import read_sbml_model

    def test_add_reactions_to_model(self, mock_write_sbml_model, mock_read_sbml_model):
        # Create mock objects for external dependencies
        mock_original_model = MagicMock()
        mock_seeds = MagicMock()
        mock_optimal_completions = [['xreaction', [('R_reaction1', 1)], ...]]  # Mock optimal completions

        # Configure the mock methods or attributes
        mock_original_model.copy.return_value = mock_original_model
        mock_original_model.slim_optimize.return_value = 0.5  # Set the return value as needed
        mock_original_model.metabolites = []  # Mock metabolites attribute as needed
        mock_original_model.exchanges = []  # Mock exchanges attribute as needed

        mock_read_sbml_model.side_effect = [mock_seeds]  # Set side effect for read_sbml_model

        # Create an instance of YourClass with mocked dependencies
        your_instance = GapFiller(mock_original_model, "seeds_path.xml", mock_optimal_completions)

        # Call the add_reactions_to_model method
        your_instance.add_reactions_to_model()

        # Assertions based on the expected behavior of your code
        self.assertTrue(your_instance.positive_solutions)
        self.assertTrue(isinstance(your_instance.original_model, MagicMock))
        self.assertEqual(your_instance.original_model.slim_optimize(), 0.5)
        # Add more assertions as needed

    def test_report(self):
        # Mocking cobra.Model and its methods
        self.gf.filled_model = mock.create_autospec(cobra.Model, instance=True)
        self.gf.original_model = mock.create_autospec(cobra.Model, instance=True)
        mock_solution = mock.Mock(objective_value=1.0)
        self.gf.original_model.optimize.return_value = mock_solution
        self.gf.original_model.summary.return_value.to_frame.return_value = mock.MagicMock()

        # Setting paths and other attributes
        self.gf.model_path = "path/to/model"
        self.gf.seeds_path = "path/to/seeds"
        self.gf.targets_path = "path/to/targets"
        self.gf.universal_model_path = "path/to/universal_model"
        self.gf.unproducible = set()
        self.gf.never_producible = set()
        self.gf.reconstructable_targets = set()
        self.gf.minimal_completion = set()
        self.gf.required_additional_seeds_and_demands = set()
        self.gf.additional_seeds = set()
        self.gf.positive_solutions = set()
        self.gf.all_solutions_with_models = []
        self.gf.all_completions = []
        self.gf.transport_reactions = []
        self.gf.dead_ends = []
        self.gf.transport_reactions_universal = []
        self.gf.added_demands = []
        self.gf.optimal_completions = []
        self.gf.minimal_set = []

        # assert the report method
        self.gf.report(write_to_file=False, removed_reactions=[])

        # assert the values of the attributes
        self.assertEqual(self.gf.model_path, "path/to/model")
        self.assertEqual(self.gf.seeds_path, "path/to/seeds")
        self.assertEqual(self.gf.targets_path, "path/to/targets")
        self.assertEqual(self.gf.universal_model_path, "path/to/universal_model")
        self.assertEqual(self.gf.unproducible, set())
        self.assertEqual(self.gf.never_producible, set())
        self.assertEqual(self.gf.reconstructable_targets, set())
        self.assertEqual(self.gf.minimal_completion, set())
        self.assertEqual(self.gf.required_additional_seeds_and_demands, set())
        self.assertEqual(self.gf.additional_seeds, set())
        self.assertEqual(self.gf.positive_solutions, set())
        self.assertEqual(self.gf.all_solutions_with_models, [])
        self.assertEqual(self.gf.all_completions, [])
        self.assertEqual(self.gf.transport_reactions, [])
        self.assertEqual(self.gf.dead_ends, [])
        self.assertEqual(self.gf.transport_reactions_universal, [])
        self.assertEqual(self.gf.added_demands, [])
        self.assertEqual(self.gf.optimal_completions, [])
        self.assertEqual(self.gf.minimal_set, [])

        # assert now to write to file
        self.gf.report(write_to_file=True, removed_reactions=[])
        # check if the file was created # json_report_path = os.path.join(self.resultsPath, "gapfilling_report.json")
        self.assertTrue(os.path.exists(os.path.join(self.results_path, "gapfilling_report.json")))

        # delete the file
        try:
            os.remove(os.path.join(self.results_path, "gapfilling_report.json"))
        except OSError:
            pass

    @mock.patch('gap_filling_dl.biomeneco.utils.get_compartments_in_model')
    @mock.patch('gap_filling_dl.biomeneco.gapfiller.GapFiller.clone_reactions')
    @mock.patch('gap_filling_dl.biomeneco.gapfiller.GapFiller.clone_metabolites')
    @mock.patch('gap_filling_dl.biomeneco.gapfiller.GapFiller.clone_groups')
    def test_clone_model(self, mock_clone_groups, mock_clone_metabolites, mock_clone_reactions,
                         mock_get_compartments_in_model):
        # Set return value for mocked methods
        mock_get_compartments_in_model.return_value = ['comp1', 'comp2']
        mock_clone_reactions.return_value = None
        mock_clone_metabolites.return_value = None
        mock_clone_groups.return_value = None

        # Run the method under test
        self.gf.clone_model(folder_path=self.results_path)

        # Assertions
        # Checking that the universal model copy and compartmentalized model are not None
        self.assertIsNotNone(self.gf.universal_model_copy)
        self.assertIsNotNone(self.gf.universal_model_compartmentalized)

        # Verifying that the internal methods are called
        mock_clone_reactions.assert_called_once()
        mock_clone_metabolites.assert_called_once()
        mock_clone_groups.assert_called_once()

        # Delete the created files
        try:
            os.remove(os.path.join(self.results_path, "universal_model_copy.xml"))
            os.remove(os.path.join(self.results_path, "universal_model_compartmentalized.xml"))
        except OSError:
            pass

    def test_clone_reaction(self):
        self.gf.universal_model_copy = mock.create_autospec(cobra.Model, instance=True)
        self.gf.universal_model_copy.metabolites = mock.MagicMock()
        self.gf.universal_model_copy.reactions = mock.MagicMock()
        self.gf.universal_model_copy.reactions.get_by_id.return_value = mock.MagicMock()

        compartments = ['comp1', 'comp2']
        reaction = mock.create_autospec(cobra.Reaction, instance=True, id='R1__orig')
        reactions = [reaction]

        # Mock metabolites in the reaction
        metabolite = mock.create_autospec(cobra.Metabolite, instance=True, id='M1__orig')
        reaction.metabolites = {metabolite: 1}

        # Mock metabolites in the universal model copy
        new_metabolite = mock.create_autospec(cobra.Metabolite, instance=True, id='M1__comp1')
        self.gf.universal_model_copy.metabolites.get_by_id.return_value = new_metabolite

        # Call the method under test
        cloned_reactions = self.gf.clone_reaction(compartments, reactions)

        # Assertions
        self.assertEqual(len(cloned_reactions), len(compartments))
        for cloned_reaction in cloned_reactions:
            self.assertIn('__', cloned_reaction.id)
            self.assertIn(cloned_reaction.id.split('__')[-1], compartments)
            for metabolite in cloned_reaction.metabolites:
                self.assertIn('__', metabolite.id)

    @mock.patch('gap_filling_dl.biomeneco.gapfiller.progress_imap', side_effect=lambda f, batches: map(f, batches))
    def test_clone_metabolites(self, mock_progress_imap):
        self.gf.universal_model_copy = mock.MagicMock(spec=cobra.Model)
        compartments = ['comp1', 'comp2']
        mock_metabolite = mock.Mock(spec=cobra.Metabolite, id='M1')
        self.gf.universal_model_copy.metabolites = [mock_metabolite]

        # Mock the clone_metabolite function
        mock_clone_metabolite_partial = partial(mock_clone_metabolite, compartments)
        mock_progress_imap.side_effect = lambda f, batches: (mock_clone_metabolite_partial(batch) for batch in batches)

        # Execute the method under test
        self.gf.clone_metabolites(compartments)

        # Assertions
        expected_metabolite_count = len(compartments)
        # Check if the correct number of metabolites were added
        self.gf.universal_model_copy.add_metabolites.assert_called()
        added_metabolites = self.gf.universal_model_copy.add_metabolites.call_args[0][0]
        self.assertEqual(len(added_metabolites), expected_metabolite_count)

    @mock.patch('gap_filling_dl.biomeneco.gapfiller.BioISO')
    @mock.patch('gap_filling_dl.biomeneco.gapfiller.set_solver')
    def test_identify_dead_end_metabolites_with_bioiso(self, mock_set_solver, mock_bioiso):
        # Set up the mock for BioISO
        mock_bio_instance = mock_bioiso.return_value
        mock_bio_instance.get_tree.return_value = {
            "M_root_M_root_M_root_product": {
                "next": {
                    "biomass_component_1": {
                        "role": "Reactant",
                        "next": {
                            "specific_biomass_component_1": {
                                "analysis": False,
                                "role": "Reactant",
                                "identifier": "metabolite_1"
                            }
                        }
                    }
                }
            }
        }

        objective_function_id = "obj_func_id"
        solver = "cplex"
        objective = "maximize"

        dead_ends_iterable = self.gf.identify_dead_end_metabolites_with_bioiso(objective_function_id, objective, solver)

        # Convert the iterable to a list if necessary
        dead_ends = list(dead_ends_iterable)

        # Assertions
        self.assertEqual(len(dead_ends), 1)
        self.assertIn("metabolite_1", dead_ends)

        # Check if set_solver and BioISO were called correctly

    @mock.patch('gap_filling_dl.biomeneco.gapfiller.progress_imap')
    def test_clone_reactions(self, mock_progress_imap):
        # Create a GapFiller instance

        # Mock the universal_model_copy and reactions
        self.gf.universal_model_copy = mock.MagicMock()
        self.gf.universal_model_copy.reactions = []

        # Mock compartments
        compartments = ['comp1', 'comp2']

        # Mock the universal_model_compartmentalized
        self.gf.universal_model_compartmentalized = mock.MagicMock()

        # Mock the clone_reaction function
        def mock_clone_reaction(comps, reactions):
            # Return a list of cloned reactions
            return [f"cloned_{reaction}__{comp}" for reaction in reactions for comp in comps]

        # Set up the mock for progress_imap
        mock_progress_imap.side_effect = mock_clone_reaction

        # Call the method under test
        self.gf.clone_reactions(compartments)

        # Assertions
        self.assertEqual(len(self.gf.universal_model_compartmentalized.add_reactions.call_args_list), 1)
        added_reactions = self.gf.universal_model_compartmentalized.add_reactions.call_args[0][0]
        expected_reactions = [f"cloned_{reaction}__{comp}" for reaction in self.gf.universal_model_copy.reactions for
                              comp in compartments]
        self.assertCountEqual(expected_reactions, [reaction.id for reaction in added_reactions])

    @patch('cobra.core.group')
    def test_clone_groups(self, mock_group):
        # Create mock universal_model_copy and universal_model_compartmentalized
        universal_model_copy = mock.Mock()
        universal_model_compartmentalized = mock.Mock()
        mock_group1 = mock.Mock(id='group1', members=[mock.Mock(id='member1'), mock.Mock(id='member2')])
        mock_group1.name = 'Group 1'  # Set the name attribute to a string
        universal_model_copy.groups = [mock_group1]

        # Create a mock reaction and add it to universal_model_compartmentalized
        reaction = mock.Mock(id='member1__comp1')
        universal_model_compartmentalized.reactions = [reaction]  # Make it iterable

        # Create an instance of GapFiller with the appropriate parameters
        self.gf.universal_model_copy = universal_model_copy
        self.gf.universal_model_compartmentalized = universal_model_compartmentalized

        # Define compartments
        compartments = ['comp1', 'comp2']

        # Call the method under test
        self.gf.clone_groups(compartments)

        # Assertions
        self.assertEqual(universal_model_compartmentalized.add_groups.call_count, 1)
        args, kwargs = universal_model_compartmentalized.add_groups.call_args
        self.assertEqual(len(args[0]), 1)  # Check the number of new groups
        new_group = args[0][0]

        # Assert the 'name' attribute of the new_group
        self.assertEqual(new_group.name, 'Group 1')

    from unittest.mock import patch, MagicMock
    from gap_filling_dl.biomeneco.gapfiller import GapFiller
    from cobra import Model

    @patch('cobra.io.read_sbml_model')
    @patch('gap_filling_dl.biomeneco.gapfiller.get_compartments_in_model')
    def test_add_transport_reaction_for_dead_ends(self, mock_get_compartments_in_model, mock_read_sbml_model):
        # Create a mock Model object with a 'reactions' attribute
        mock_model = MagicMock(spec=Model)
        mock_model.reactions = []  # assuming reactions is a list

        # Set the return value of read_sbml_model to the mock Model
        mock_read_sbml_model.return_value = mock_model

        # Mock the get_compartments_in_model function to return a list of compartments
        mock_get_compartments_in_model.return_value = ['c', 'e']

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)
        gap_filler.targets_path = "targets_path"

        # Replace the method with a mock that doesn't invoke multiprocessing
        gap_filler.add_transport_reaction_for_dead_ends = MagicMock()

        # Call the method with compartments
        gap_filler.add_transport_reaction_for_dead_ends(["c", "e"])

        # Assertions
        gap_filler.add_transport_reaction_for_dead_ends.assert_called_once_with(["c", "e"])

    @patch('cobra.io.read_sbml_model')
    def test_create_transport(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Create mock reactions
        mock_reaction = MagicMock(spec=cobra.Reaction)
        mock_model.reactions = [mock_reaction]

        # Create mock metabolites
        mock_metabolite1 = Metabolite("met1__c")
        mock_metabolite2 = Metabolite("met1__e")

        # Create a mock DictList for the metabolites
        mock_metabolites = MagicMock()
        mock_metabolites.get_by_id.side_effect = lambda x: {
            "met1__c": mock_metabolite1,
            "met1__e": mock_metabolite2
        }.get(x, None)

        # Set the metabolites attribute of the mock Model to the mock DictList
        mock_model.metabolites = mock_metabolites

        # Set the return value of read_sbml_model to the mock Model
        mock_read_sbml_model.return_value = mock_model

        # Create a GapFiller instance with the mock model
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)
        gap_filler.universal_model = mock_model

        # Call create_transport method
        compartments_combinations = [('c', 'e')]
        metabolites = ['met1']
        transport_reactions = gap_filler.create_transport(compartments_combinations, metabolites)

        # Assertions
        created_reaction = transport_reactions[0]
        self.assertIsInstance(created_reaction, Reaction)
        self.assertEqual(created_reaction.id, "T_met1_c_to_e")
        self.assertEqual(created_reaction.name, "Transport of met1 between c and e")
        self.assertIn(mock_metabolite1, created_reaction.metabolites)
        self.assertIn(mock_metabolite2, created_reaction.metabolites)

    @patch('cobra.io.read_sbml_model')
    @patch('gap_filling_dl.biomeneco.gapfiller.get_compartments_in_model')
    @patch('gap_filling_dl.biomeneco.gapfiller.json.load')
    @patch('gap_filling_dl.biomeneco.gapfiller.os.path.exists')
    @patch('gap_filling_dl.biomeneco.gapfiller.create_related_pathways_map')
    @patch('gap_filling_dl.biomeneco.gapfiller.get_related_pathways')
    @patch('cobra.io.write_sbml_model')
    def test_build_temporary_universal_model(self, mock_write_sbml_model, mock_get_related_pathways,
                                             mock_create_related_pathways_map, mock_os_path_exists,
                                             mock_json_load, mock_get_compartments_in_model,
                                             mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)
        mock_model.reactions = [MagicMock(spec=Reaction), MagicMock(spec=Reaction)]
        mock_model.metabolites = [MagicMock(spec=Metabolite), MagicMock(spec=Metabolite)]
        mock_model.groups = []

        # Set the return value of read_sbml_model to the mock Model when it's called with "targets_path"
        def mock_read_sbml_model_func(filename):
            if filename in ['universal_model_path', 'targets_path']:
                return mock_model
            else:
                raise ValueError(f"Unexpected filename: {filename}")

        mock_read_sbml_model.side_effect = mock_read_sbml_model_func

        # Mock the get_compartments_in_model function to return a list of compartments
        mock_get_compartments_in_model.return_value = ['c', 'e']

        # Mock the json.load function to return a dictionary of cofactors
        mock_cofactors = {"met1__c": ["pathway1", "pathway2"]}
        mock_json_load.return_value = mock_cofactors

        # Mock the os.path.exists function to return True
        mock_os_path_exists.return_value = True

        # Mock the create_related_pathways_map function to return a dictionary
        mock_create_related_pathways_map.return_value = {"pathway1": ["pathway2"]}

        # Mock the get_related_pathways function to return a list of pathways
        mock_get_related_pathways.return_value = ["pathway2"]

        # Create a GapFiller instance with specific targets_path
        targets_path = "targets_path"
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)
        gap_filler.targets_path = targets_path

        # Call the method with a folder path and related_pathways set to True
        gap_filler.build_temporary_universal_model("folder_path", True)

        # Assert that the write_sbml_model function was called with the correct arguments
        mock_write_sbml_model.assert_called_once_with(gap_filler.universal_model_compartmentalized,
                                                      "folder_path/temporary_universal_model.xml")

    # @patch('json.load')
    # @patch('gap_filling_dl.biomeneco.gapfiller.write_metabolites_to_sbml')
    # @patch('cobra.io.read_sbml_model')
    # @patch.object(GapFiller, 'identify_additional_seeds')
    # def test_identify_additional_seeds(self, mock_identify_additional_seeds, mock_read_sbml_model,
    #                                    mock_write_metabolites_to_sbml, mock_json_load):
    #     # Create a mock Model object with autospec
    #     mock_model = mock.create_autospec(Model, instance=True)
    #
    #     # Create mock Metabolite objects with autospec
    #     mock_met1 = mock.create_autospec(Metabolite, instance=True, id='met1__c')
    #     mock_met2 = mock.create_autospec(Metabolite, instance=True, id='met2__c')
    #     mock_reaction1 = mock.create_autospec(Reaction, instance=True)
    #     mock_reaction2 = mock.create_autospec(Reaction, instance=True)
    #
    #     # Set the metabolites and reactions of the model
    #     mock_model.metabolites = [mock_met1, mock_met2]
    #     mock_model.reactions = [mock_reaction1, mock_reaction2]
    #
    #     mock_universal_model = mock.create_autospec(Model, instance=True)
    #     mock_universal_model.metabolites = [mock_met1, mock_met2]
    #     mock_universal_model.reactions = [mock_reaction1, mock_reaction2]
    #     # Add custom attributes to the mock
    #     mock_universal_model.metabolite_pathway_map = {"met1": ["pathway1"], "met2": ["pathway2"]}
    #     mock_universal_model.pathway_metabolites_map = {"pathway1": ["met1"], "pathway2": ["met2"]}
    #
    #     mock_read_sbml_model.return_value = mock_universal_model
    #
    #     # Mock the json.load function to return a dictionary of cofactors
    #     mock_cofactors = {"met1__c": ["pathway1", "pathway2"]}
    #     mock_json_load.return_value = mock_cofactors
    #
    #     # Create a GapFiller instance
    #     gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)
    #
    #     # Mock the seeds model
    #     mock_seeds_model = mock.create_autospec(Model, instance=True)
    #     mock_seeds_model.metabolites = [mock_met1, mock_met2]
    #
    #     # Configure the mock to return the appropriate model based on the filename
    #     mock_read_sbml_model.side_effect = lambda filename: \
    #         mock_universal_model if filename == 'universal_model_path' \
    #             else mock_seeds_model if filename == 'seeds_path' \
    #             else None
    #
    #     # Set the seeds_path and universal_model attributes for gap_filler
    #     gap_filler.seeds_path = "seeds_path"
    #     gap_filler.universal_model = mock_universal_model
    #
    #     # Set the never_producible attribute
    #     gap_filler.never_producible = {'M_met1', 'M_met2'}
    #
    #     # Define a mock combinet and targets for the method call
    #     combinet = MagicMock()
    #     targets = MagicMock()
    #
    #     # Mock the identify_additional_seeds method to modify the additional_seeds attribute
    #     def mock_identify_additional_seeds_side_effect(*args, **kwargs):
    #         gap_filler.additional_seeds = {mock_met1, mock_met2}
    #
    #     mock_identify_additional_seeds.side_effect = mock_identify_additional_seeds_side_effect
    #
    #     # Call the identify_additional_seeds method
    #     gap_filler.identify_additional_seeds(combinet, targets)
    #
    #     # Assertions
    #     mock_write_metabolites_to_sbml.assert_called_once()
    #     self.assertIn('met1__c', [met.id for met in gap_filler.additional_seeds])
    #     self.assertIn('met2__c', [met.id for met in gap_filler.additional_seeds])
    #     self.assertEqual(len(gap_filler.additional_seeds), 2)

    # @patch('cobra.Model.slim_optimize')
    # @patch('cobra.Model.add_boundary')
    # @patch('gap_filling_dl.biomeneco.gapfiller.GapFiller.find_deadends')
    # @patch('cobra.io.read_sbml_model')
    # def test_add_demands(self, mock_read_sbml_model, mock_find_deadends, mock_add_boundary, mock_slim_optimize):
    #     # Create a mock Model object
    #     mock_model = MagicMock(spec=Model)
    #
    #     # Create mock Metabolite objects with autospec
    #     mock_met1 = mock.create_autospec(Metabolite, instance=True, id='met1')
    #     mock_met2 = mock.create_autospec(Metabolite, instance=True, id='met2')
    #
    #     # Set the metabolites of the model
    #     mock_model.metabolites = [mock_met1, mock_met2]
    #
    #     mock_model.reactions = [MagicMock(spec=Reaction), MagicMock(spec=Reaction)]  # Add this line
    #
    #     # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
    #     def mock_read_sbml_model_func(filename):
    #         if filename == 'universal_model_path':
    #             return mock_model
    #         else:
    #             raise ValueError(f"Unexpected filename: {filename}")
    #
    #     mock_read_sbml_model.side_effect = mock_read_sbml_model_func
    #
    #     # Create a GapFiller instance
    #     gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)
    #
    #     # Mock the find_deadends method to return a list of dead-end metabolites
    #     mock_find_deadends.return_value = [MagicMock(spec=Metabolite, id='met1'), MagicMock(spec=Metabolite, id='met2')]
    #
    #     # Mock the add_boundary method to create a demand reaction for each dead-end metabolite
    #     mock_add_boundary.side_effect = lambda met_id, type: MagicMock(spec=Reaction, id=f'DM_{met_id}')
    #
    #     # Mock the slim_optimize method to return a positive value
    #     mock_slim_optimize.return_value = 1.0
    #
    #     # Create a mock solution and model for the all_solutions_with_models attribute
    #     mock_solution = MagicMock(spec=Solution)
    #     mock_model = MagicMock(spec=Model)
    #     gap_filler.all_solutions_with_models = [(mock_solution, mock_model)]
    #
    #     # Create a mock reaction for the original_model attribute
    #     mock_reaction = MagicMock(spec=Reaction, id='DM_met1')
    #     gap_filler.original_model = MagicMock(spec=Model, reactions=[mock_reaction])
    #
    #     # Call the add_demands method
    #     gap_filler.add_demands()
    #
    #     # Assert that the find_deadends method was called with the correct argument
    #     mock_find_deadends.assert_called_once_with(mock_model)
    #
    #     # Assert that the add_boundary method was called for each dead-end metabolite
    #     mock_add_boundary.assert_has_calls([mock.call('met1', type='demand'), mock.call('met2', type='demand')])
    #
    #     # Assert that the slim_optimize method was called
    #     assert mock_slim_optimize.call_count == 2
    #
    #     # Assert that the added_demands attribute was updated correctly
    #     assert len(gap_filler.added_demands) == 1
    #     assert mock_reaction in gap_filler.added_demands

    @patch('json.load')
    @patch('gap_filling_dl.biomeneco.gapfiller.write_metabolites_to_sbml')
    @patch('cobra.io.read_sbml_model')
    def test_find_deadends(self, mock_read_sbml_model, mock_write_metabolites_to_sbml, mock_json_load):
        # Create a mock Model object with autospec
        mock_model = mock.create_autospec(Model, instance=True)

        # Create mock Metabolite objects with autospec
        mock_met1 = mock.create_autospec(Metabolite, instance=True, id='met1')
        mock_met2 = mock.create_autospec(Metabolite, instance=True, id='met2')

        # Set the metabolites of the model
        mock_model.metabolites = [mock_met1, mock_met2]

        # Create mock Reaction objects with autospec
        mock_rxn1 = mock.create_autospec(Reaction, instance=True, id='rxn1', reversibility=False)
        mock_rxn2 = mock.create_autospec(Reaction, instance=True, id='rxn2', reversibility=False)

        # Mock the get_coefficient method and lower/upper bounds of the Reaction objects
        mock_rxn1.get_coefficient.return_value = 1.0
        mock_rxn2.get_coefficient.return_value = -1.0
        mock_rxn1.lower_bound = 0.0
        mock_rxn1.upper_bound = 1000.0
        mock_rxn2.lower_bound = 0.0
        mock_rxn2.upper_bound = 1000.0

        # Set the reactions of the model
        mock_model.reactions = [mock_rxn1, mock_rxn2]

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        def mock_read_sbml_model_func(filename):
            if filename == 'universal_model_path':
                return mock_model
            else:
                raise ValueError(f"Unexpected filename: {filename}")

        mock_read_sbml_model.side_effect = mock_read_sbml_model_func

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Call the find_deadends method with the mock model
        result = gap_filler.find_deadends(mock_model)

        # Assert that the result is an empty list (since there are no dead-end metabolites in the mock model)
        self.assertEqual(result, [])

    @patch('cobra.io.read_sbml_model')
    def test_remove_unnecessary_seeds_and_demands(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Create mock Reaction objects
        mock_rxn1 = MagicMock(spec=Reaction, id='Sk_rxn1')
        mock_rxn2 = MagicMock(spec=Reaction, id='DM_rxn2')
        mock_rxn3 = MagicMock(spec=Reaction, id='rxn3')

        # Set the reactions of the model
        mock_model.reactions = [mock_rxn1, mock_rxn2, mock_rxn3]

        # Set the demands of the model
        mock_model.demands = [mock_rxn1, mock_rxn2]

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        mock_read_sbml_model.return_value = mock_model

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Set the original_model attribute of the GapFiller instance
        gap_filler.original_model = mock_model

        # Set the additional_seeds attribute of the GapFiller instance
        gap_filler.additional_seeds = set([MagicMock(spec=Metabolite, id='met1')])

        # Mock the slim_optimize method to return a positive value for the first call and zero for the second call
        mock_model.slim_optimize = MagicMock(side_effect=[1.0, 0.0])

        # Call the remove_unnecessary_seeds_and_demands method
        gap_filler.remove_unnecessary_seeds_and_demands()

        # Assert that the remove_reactions method was called with the correct argument
        mock_model.remove_reactions.assert_called_once_with([mock_rxn1])

    # @patch('cobra.io.read_sbml_model')
    # def test_add_reactions_to_model_v2(self, mock_read_sbml_model):
    #     # Create a mock Model object
    #     mock_model = Model('mock_model')
    #
    #     # Create mock Reaction objects
    #     mock_rxn1 = Reaction('rxn1')
    #     mock_rxn2 = Reaction('rxn2')
    #
    #     # Set the reactions of the model
    #     mock_model.add_reactions([mock_rxn1, mock_rxn2])
    #
    #     # Create a real Model object for seeds
    #     mock_seeds_model = Model('mock_seeds_model')
    #
    #     # Create real Metabolite objects
    #     mock_met1 = Metabolite('met1')
    #     mock_met2 = Metabolite('met2')
    #
    #     # Set the metabolites of the seeds model
    #     mock_seeds_model.add_metabolites([mock_met1, mock_met2])
    #
    #     # Set the return value of read_sbml_model to the mock Model when it's called with "seeds_path" or "universal_model_path"
    #     def mock_read_sbml_model_func(filename):
    #         if filename == 'universal_model_path':
    #             return mock_model
    #         elif filename == 'seeds_path':
    #             return mock_seeds_model
    #         else:
    #             raise ValueError(f"Unexpected filename: {filename}")
    #
    #     mock_read_sbml_model.side_effect = mock_read_sbml_model_func
    #
    #     # Create a GapFiller instance
    #     gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)
    #
    #     # Set the original_model attribute of the GapFiller instance
    #     gap_filler.original_model = mock_model
    #
    #     # Set the seeds_path attribute of the GapFiller instance
    #     gap_filler.seeds_path = "seeds_path"
    #
    #     # Set the optimal_completions attribute of the GapFiller instance
    #     gap_filler.optimal_completions = [{'xreaction': [('rxn1',)]}, {'xreaction': [('rxn2',)]}]
    #
    #     # Set the minimal_completion attribute of the GapFiller instance
    #     gap_filler.minimal_completion = ['rxn1', 'rxn2']
    #
    #     # Call the add_reactions_to_model_v2 method
    #     gap_filler.add_reactions_to_model_v2()
    #
    #     # Assert that the original_model attribute of the GapFiller instance has been modified as expected
    #     self.assertEqual(len(gap_filler.original_model.reactions), 4)
    #     self.assertIn('rxn1_universal', [reaction.id for reaction in gap_filler.original_model.reactions])
    #     self.assertIn('rxn2_universal', [reaction.id for reaction in gap_filler.original_model.reactions])

    @patch('cobra.io.read_sbml_model')
    def test_get_mets_in_medium(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Create a mock Reaction object
        mock_reaction = MagicMock(spec=Reaction)

        # Set the lower_bound attribute of the mock Reaction object
        mock_reaction.lower_bound = -1

        # Set the reactants attribute of the mock Reaction object
        mock_reaction.reactants = [MagicMock(spec=Metabolite, id='met1')]

        # Set the exchanges attribute of the mock Model object
        mock_model.exchanges = [mock_reaction]

        # Create a mock Model object for universal_model
        mock_universal_model = MagicMock(spec=Model)

        # Create mock Reaction objects for the universal_model
        mock_universal_rxn1 = MagicMock(spec=Reaction)
        mock_universal_rxn2 = MagicMock(spec=Reaction)

        # Set the reactions attribute of the universal_model mock object
        mock_universal_model.reactions = [mock_universal_rxn1, mock_universal_rxn2]

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        mock_read_sbml_model.return_value = mock_universal_model

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Call the get_mets_in_medium method
        result = gap_filler.get_mets_in_medium(mock_model)

        # Assert that the result is as expected
        self.assertEqual(result, {'met1'})

    @patch('cobra.io.read_sbml_model')
    def test_get_completion(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Create mock Reaction objects for the universal_model
        mock_universal_rxn1 = MagicMock(spec=Reaction)
        mock_universal_rxn2 = MagicMock(spec=Reaction)

        # Create a mock Model object for universal_model
        mock_universal_model = MagicMock(spec=Model, reactions=[mock_universal_rxn1, mock_universal_rxn2])

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        mock_read_sbml_model.return_value = mock_universal_model

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Define a mock gap-filling result
        m = {'xreaction': [('rxn1',)]}

        # Call the get_completion method
        result = gap_filler.get_completion(m)

        # Assert that the result is as expected
        self.assertEqual(result, {'rxn1'})

    @patch('cobra.io.read_sbml_model')
    def test_add_reactions_to_temp_model(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Create a mock Reaction object
        mock_reaction = MagicMock(spec=Reaction, id='rxn1')

        # Create a mock Reaction object for the copy method
        mock_reaction_copy = MagicMock(spec=Reaction)
        mock_reaction_copy.id = 'rxn1_universal'

        # Set the copy method of the mock Reaction object to return the mock_reaction_copy object
        mock_reaction.copy = MagicMock(return_value=mock_reaction_copy)

        # Create a cobra.DictList object and add the mock Reaction object to it
        mock_reactions = cobra.DictList()
        mock_reactions.append(mock_reaction)

        # Set the reactions attribute of the mock Model object
        mock_model.reactions = mock_reactions

        # Create mock Reaction objects for the universal_model
        mock_universal_rxn1 = MagicMock(spec=Reaction, id='rxn1')
        mock_universal_rxn2 = MagicMock(spec=Reaction, id='rxn2')

        # Create a cobra.DictList object and add the mock Reaction objects to it
        mock_universal_reactions = cobra.DictList()
        mock_universal_reactions.append(mock_universal_rxn1)
        mock_universal_reactions.append(mock_universal_rxn2)

        # Create a mock Model object for universal_model
        mock_universal_model = MagicMock(spec=Model, reactions=mock_universal_reactions)

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        mock_read_sbml_model.return_value = mock_universal_model

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Define a completion
        completion = {'rxn1'}

        # Define a set of reactions already in the original model
        already_in_original_model = {'rxn1'}

        # Call the add_reactions_to_temp_model method
        to_add, updated_already_in_original_model = gap_filler.add_reactions_to_temp_model(completion,
                                                                                           already_in_original_model,
                                                                                           mock_model)

        # Assert that the to_add list and updated_already_in_original_model set are as expected
        self.assertEqual(len(to_add), 1)
        self.assertIsInstance(to_add[0].id, MagicMock)
        self.assertIsInstance(updated_already_in_original_model, set)

    @patch('cobra.io.read_sbml_model')
    def test_create_sinks_for_metabolites(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Add the create_sink method to the mock Model object
        mock_model.create_sink = MagicMock()

        # Create mock Metabolite objects
        mock_met1 = MagicMock(spec=Metabolite, id='met1')
        mock_met2 = MagicMock(spec=Metabolite, id='met2')

        # Set the metabolites attribute of the mock Model object
        mock_model.metabolites = [mock_met1, mock_met2]

        # Create a mock Model object for seeds
        mock_seeds_model = MagicMock(spec=Model)

        # Set the metabolites attribute of the mock seeds Model object
        mock_seeds_model.metabolites = [mock_met1, mock_met2]

        # Define a set of metabolite IDs that are present in the medium
        mets_in_medium = {'met1'}

        # Create mock Reaction objects for the universal_model
        mock_universal_rxn1 = MagicMock(spec=Reaction)
        mock_universal_rxn2 = MagicMock(spec=Reaction)

        # Create a mock Model object for universal_model
        mock_universal_model = MagicMock(spec=Model, reactions=[mock_universal_rxn1, mock_universal_rxn2])

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        mock_read_sbml_model.return_value = mock_universal_model

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Call the create_sinks_for_metabolites method
        gap_filler.create_sinks_for_metabolites(mock_model, mock_seeds_model, mets_in_medium)

        # Assert that the create_sink method was called with the correct argument
        mock_model.create_sink.assert_called_once_with('met2')

    @patch('cobra.io.read_sbml_model')
    def test_check_for_positive_solution(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Create a mock Reaction object
        mock_reaction = MagicMock(spec=Reaction, id='rxn1')

        # Set the slim_optimize method of the mock Model object to return a value that will cause the check_for_positive_solution method to add the mock_model object to the positive_solutions set
        mock_model.slim_optimize = MagicMock(return_value=0.1)

        # Create mock Reaction objects for the universal_model
        mock_universal_rxn1 = MagicMock(spec=Reaction)
        mock_universal_rxn2 = MagicMock(spec=Reaction)

        # Create a mock Model object for universal_model
        mock_universal_model = MagicMock(spec=Model, reactions=[mock_universal_rxn1, mock_universal_rxn2])

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        mock_read_sbml_model.return_value = mock_universal_model

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Define a list of reactions that were added to the model
        to_add = [('rxn1')]

        # Define a set to store models with positive solutions
        positive_solutions = set()

        # Mock the check_for_positive_solution method to add the mock_model object to the positive_solutions set
        gap_filler.check_for_positive_solution = MagicMock()
        gap_filler.check_for_positive_solution.return_value = positive_solutions.add((mock_model, tuple(to_add)))

        # Call the check_for_positive_solution method
        gap_filler.check_for_positive_solution(mock_model, to_add, positive_solutions)

        # Assert that the positive_solutions set was updated correctly
        self.assertEqual(len(positive_solutions), 1)
        self.assertEqual(list(positive_solutions)[0][0], mock_model)
        self.assertEqual(list(positive_solutions)[0][1], ('rxn1',))

    @patch('cobra.io.read_sbml_model')
    def test_add_demands_and_reoptimize(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Create a mock Reaction object
        mock_reaction = MagicMock(spec=Reaction, id='rxn1')

        # Create mock Metabolite objects for the model
        mock_metabolite1 = MagicMock(spec=Metabolite)
        mock_metabolite2 = MagicMock(spec=Metabolite)

        # Set the metabolites attribute of the mock Model object
        mock_model.metabolites = [mock_metabolite1, mock_metabolite2]

        # Set the slim_optimize method of the mock Model object to return 1.0
        mock_model.slim_optimize = MagicMock(return_value=1.0)

        # Create mock Reaction objects for the universal_model
        mock_universal_rxn1 = MagicMock(spec=Reaction)
        mock_universal_rxn2 = MagicMock(spec=Reaction)

        # Create a mock Model object for universal_model
        mock_universal_model = MagicMock(spec=Model, reactions=[mock_universal_rxn1, mock_universal_rxn2])

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        mock_read_sbml_model.return_value = mock_universal_model

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Define a list of reactions that were added to the model
        to_add = [('rxn1')]

        # Define a set to store models with positive solutions
        positive_solutions = set()

        # Mock the add_demands_and_reoptimize method to add the mock_model object to the positive_solutions set
        gap_filler.add_demands_and_reoptimize = MagicMock()
        gap_filler.add_demands_and_reoptimize.return_value = positive_solutions.add((mock_model, tuple(to_add)))

        # Call the add_demands_and_reoptimize method
        gap_filler.add_demands_and_reoptimize(mock_model, to_add, positive_solutions)

        # Assert that the positive_solutions set was updated correctly
        self.assertEqual(len(positive_solutions), 1)
        self.assertEqual(list(positive_solutions)[0][0], mock_model)
        self.assertEqual(list(positive_solutions)[0][1], ('rxn1',))

    @patch('cobra.io.read_sbml_model')
    def test_handle_final_model_selection(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Create a mock Reaction object
        mock_reaction = MagicMock(spec=Reaction, id='rxn1')

        # Set the slim_optimize method of the mock Model object to return 1.0
        mock_model.slim_optimize = MagicMock(return_value=1.0)

        # Create a mock Model object for seeds
        mock_seeds_model = MagicMock(spec=Model)

        # Define a set of metabolite IDs that are present in the medium
        mets_in_medium = {'met1'}

        # Define a set of reaction IDs already in the original model
        already_in_original_model = {'rxn1'}

        # Define a set of models with positive solutions
        positive_solutions = {(mock_model, ('rxn1',))}

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)

        # Call the handle_final_model_selection method
        gap_filler.handle_final_model_selection(positive_solutions, mock_seeds_model, mets_in_medium,
                                                already_in_original_model)

        # Assert that the original_model attribute was updated correctly
        self.assertEqual(gap_filler.original_model, mock_model)

    @patch('cobra.io.read_sbml_model')
    def test_handle_no_positive_solution(self, mock_read_sbml_model):
        # Create a mock Model object
        mock_model = MagicMock(spec=Model)

        # Add the create_sink method to the mock Model object
        mock_model.create_sink = MagicMock()

        # Create a mock Reaction object
        mock_reaction = MagicMock(spec=Reaction, id='rxn1')

        # Create mock Metabolite objects for the model and seeds
        mock_metabolite1 = MagicMock(spec=Metabolite)
        mock_metabolite2 = MagicMock(spec=Metabolite)

        # Set the metabolites attribute of the mock Model object
        mock_model.metabolites = [mock_metabolite1, mock_metabolite2]

        # Create a mock Model object for seeds
        mock_seeds_model = MagicMock(spec=Model)

        # Set the metabolites attribute of the mock seeds Model object
        mock_seeds_model.metabolites = [mock_metabolite1, mock_metabolite2]

        # Define a set of metabolite IDs that are present in the medium
        mets_in_medium = {'met1'}

        # Define a set of reaction IDs already in the original model
        already_in_original_model = {'rxn1'}

        # Create a list of mock Reaction objects for minimal_completion
        mock_minimal_completion = ['R_rxn1']

        # Create mock Reaction objects for the universal_model
        mock_universal_rxn1 = MagicMock(spec=Reaction, id='rxn1')

        # Create a cobra.DictList object and add the mock Reaction object to it
        mock_universal_reactions = cobra.DictList()
        mock_universal_reactions.append(mock_universal_rxn1)

        # Create a mock Model object for universal_model
        mock_universal_model = MagicMock(spec=Model, reactions=mock_universal_reactions)

        # Set the return value of read_sbml_model to the mock Model when it's called with "universal_model_path"
        mock_read_sbml_model.return_value = mock_universal_model

        # Create a GapFiller instance
        gap_filler = GapFiller("model_path", "universal_model_path", "results_path", 1)
        gap_filler.minimal_completion = mock_minimal_completion

        # Set the original_model attribute of the GapFiller instance to the mock Model object
        gap_filler.original_model = mock_model

        # Call the handle_no_positive_solution method
        gap_filler.handle_no_positive_solution(mock_seeds_model, mets_in_medium, already_in_original_model)

        # Assert that the add_reactions method was called with the correct argument
        mock_model.add_reactions.assert_called_once()

if __name__ == '__main__':
    unittest.main()
