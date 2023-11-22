import glob
import shutil
import tempfile
import unittest
import os

from cobra import DictList
import cobra.io
from cobra import Reaction, Metabolite
from cobra.core import Group
from mock.mock import create_autospec, MagicMock, PropertyMock
from unittest.mock import patch

from pandas import DataFrame

from src.gap_filling_dl.biomeneco.model import Model
from bioiso import load
from tests import TEST_DIR


class TestModel(unittest.TestCase):
    def setUp(self):
        # Create a mock cobra model for testing
        self.cobra_model = MagicMock(spec=cobra.Model)
        self.cobra_model.id = "test_model"

        # Create mock reactions, metabolites, and groups
        mock_reaction = create_autospec(Reaction, instance=True)
        mock_metabolite = create_autospec(Metabolite, instance=True)
        mock_group = create_autospec(Group, instance=True)

        # Add the mock reactions, metabolites, and groups to the mock model
        self.cobra_model.reactions = [mock_reaction]
        self.cobra_model.metabolites = [mock_metabolite]
        self.cobra_model.groups = [mock_group]

        # Add a mock genes attribute to the mock model
        self.cobra_model.genes = []

        # mock the id attribute of the cobra.Model instance
        type(self.cobra_model).id = PropertyMock(return_value="test_model")

        self.objective_function_id = "test_objective_function"

    def test_init(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        self.assertIsNone(model.id)
        self.assertEqual(model.objective_function_id, self.objective_function_id)

        self.assertEqual(model.objective_function_id, self.objective_function_id)
        self.assertIsNone(model.e_res_precursors)
        self.assertIsNone(model.precursors_reactions)
        self.assertIsNotNone(model.reaction_pathway_map)
        self.assertIsNotNone(model.metabolite_pathway_map)
        self.assertEqual(model.pathway_metabolites_map, {})
        self.assertIsNone(model.pre_precursors)
        self.assertIsNone(model.bio_precursors)
        self.assertIsNone(model._seeds)
        self.assertIsNone(model._targets)
        self.assertEqual(model.model_old, [])

    def test_from_sbml(self):
        # Mock the cobra.io.read_sbml_model function to return the mock cobra model
        with patch('cobra.io.read_sbml_model', return_value=self.cobra_model):
            # Call the from_sbml method
            model = Model.from_sbml('mock_path', self.objective_function_id)

        # Check that the returned Model instance has the correct attributes
        self.assertIsNone(model.id)
        self.assertEqual(model.objective_function_id, self.objective_function_id)

    def test_str(self):
        # Mock the summary method to return a specific string
        self.cobra_model.summary = MagicMock(return_value="test_summary")

        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Call the __str__ method and check the returned string
        self.assertEqual(str(model), "test_summary")

    def test_seeds(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Set the _seeds attribute of the model
        model._seeds = [('metabolite1', 'compartment1'), ('metabolite2', 'compartment2')]

        # Check that the seeds property returns the correct value
        self.assertEqual(model.seeds, [('metabolite1', 'compartment1'), ('metabolite2', 'compartment2')])

    def test_seeds_property(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Set the seeds property
        seeds = [('metabolite1', 'compartment1'), ('metabolite2', 'compartment2')]
        model.seeds = seeds

        # Check that the seeds property returns the correct value
        self.assertEqual(model.seeds, seeds)

    def test_targets_property(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Set the targets property
        targets = [('metabolite3', 'compartment3'), ('metabolite4', 'compartment4')]
        model.targets = targets

        # Check that the targets property returns the correct value
        self.assertEqual(model.targets, targets)

    def test_identify_seeds(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Mock a boundary reaction with a reactant and a lower bound less than 0
        mock_reaction = MagicMock(spec=Reaction)
        mock_reaction.lower_bound = -1
        # add the reaction to the mock 'model'
        model.reactions = [mock_reaction]

        mock_metabolite = MagicMock(spec=Metabolite)
        mock_metabolite.id = 'metabolite1'
        mock_metabolite.compartment = 'compartment1'
        mock_reaction.reactants = [mock_metabolite]
        # add the metabolite to the mock 'model'
        model.metabolites = [mock_metabolite]

        # Mock the boundary attribute to return the mock reaction when accessed
        type(self.cobra_model).boundary = PropertyMock(return_value=[mock_reaction])

        # Call the identify_seeds method
        seeds = model.identify_seeds()

        # Check that the identify_seeds method returns the correct value
        self.assertEqual(seeds, [('metabolite1', 'compartment1')])

        # Check that the seeds attribute of the Model instance is correctly set
        self.assertEqual(model.seeds, [('metabolite1', 'compartment1')])

    def test_add_seeds(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Set the seeds attribute of the Model instance
        model.seeds = [('metabolite1', 'compartment1')]

        # Define new seeds to add
        new_seeds = [('metabolite2', 'compartment2'), ('metabolite3', 'compartment3')]

        # Call the add_seeds method
        model.add_seeds(new_seeds)

        # Check that the seeds attribute of the Model instance is correctly updated
        self.assertEqual(model.seeds, [('metabolite1', 'compartment1'), ('metabolite2', 'compartment2'),
                                       ('metabolite3', 'compartment3')])

    # def test_identify_targets(self):
    #     # Initialize the Model
    #     model = Model(self.cobra_model, self.objective_function_id)
    #
    #     # Mock a reaction for the objective function
    #     mock_reaction = MagicMock(spec=Reaction)
    #     model.get_reaction = MagicMock(return_value=mock_reaction)
    #
    #     # Mock the BioISO class and its methods
    #     with patch('src.gap_filling_dl.biomeneco.model.BioISO') as mock_bioiso:
    #         instance = mock_bioiso.return_value
    #         instance.run.return_value = None
    #         instance.get_tree.return_value = {
    #             "M_root_M_root_M_root_product": {
    #                 "next": {
    #                     ('metabolite1', 'compartment1'): {
    #                         "role": "target",
    #                         "analysis": None
    #                     }
    #                 }
    #             }
    #         }
    #
    #         # Set the bio_precursors attribute to a list of real Metabolite objects
    #         metabolite = Metabolite('metabolite1', compartment='compartment1')
    #         model.bio_precursors = [metabolite]
    #
    #         # Mock the evaluate_preprecursors method to return a specific value
    #         model.evaluate_preprecursors = MagicMock(return_value=True)
    #
    #         # Call the identify_targets method
    #         targets = model.identify_targets()
    #
    #         # Check that the identify_targets method returns the correct value
    #         self.assertEqual(targets, [('metabolite1', 'compartment1')])
    #
    #         # Check that the targets attribute of the Model instance is correctly set
    #         self.assertEqual(model.targets, [('metabolite1', 'compartment1')])

    def test_validate_inputs(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Test with valid inputs
        try:
            model._validate_inputs('maximize', 'cplex')
        except ValueError:
            self.fail("_validate_inputs raised ValueError unexpectedly!")

        # Test with invalid objective
        with self.assertRaises(ValueError):
            model._validate_inputs(None, 'cplex')

        # Test with invalid solver
        with self.assertRaises(ValueError):
            model._validate_inputs('maximize', None)

        # Test with invalid objective_function_id
        model.objective_function_id = None
        with self.assertRaises(ValueError):
            model._validate_inputs('maximize', 'cplex')

    def test_extract_targets(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Define a dictionary of biomass components
        biomass_components = {
            ('metabolite1', 'compartment1'): {
                "role": "Reactant",
                "analysis": None,
                "identifier": 'metabolite1',
                "compartment": 'compartment1',
                "next": None
            },
            ('metabolite3', 'compartment3'): {
                "role": "Reactant",
                "analysis": None,
                "identifier": 'metabolite3',
                "compartment": 'compartment3',
                "next": None
            }
        }

        # Call the _extract_targets method
        targets = model._extract_targets(biomass_components)

        # Check that the _extract_targets method returns the correct value
        expected_targets = {('metabolite1', 'compartment1'), ('metabolite3', 'compartment3')}
        self.assertEqual(targets, expected_targets)

    # def test_test_e_precursors(self):
    #     # Initialize the Model
    #     model = Model(self.cobra_model, self.objective_function_id)
    #
    #     # Set the bio_precursors attribute to a list of real Metabolite objects
    #     metabolite = Metabolite('metabolite1', compartment='compartment1')
    #     model.bio_precursors = [metabolite]
    #
    #     # Mock the evaluate_preprecursors method to return a specific value
    #     model.evaluate_preprecursors = MagicMock(return_value=(0, 'metabolite1'))
    #
    #     # Call the test_e_precursors method
    #     essential_precursors = model.test_e_precursors()
    #
    #     # Check that the test_e_precursors method returns the correct value
    #     expected_essential_precursors = DataFrame({"Flux": [0]}, index=['metabolite1']).index
    #     self.assertEqual(essential_precursors, expected_essential_precursors)

    def test_evaluate_preprecursors(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Mock the evaluate_precursor method to return a specific value
        model.evaluate_precursor = MagicMock(return_value=1)

        # Set the pre_precursors attribute to a specific value
        metabolite1 = Metabolite('metabolite1', compartment='compartment1')
        metabolite2 = Metabolite('metabolite2', compartment='compartment2')
        model.pre_precursors = {'metabolite1': [metabolite2]}

        # Call the evaluate_preprecursors method
        e_precursors_res, meta_val = model.evaluate_preprecursors(['metabolite3'], metabolite1)

        # Check that the evaluate_preprecursors method returns the correct value
        self.assertEqual(e_precursors_res, {"Flux": [1]})
        self.assertEqual(meta_val, ['metabolite2'])

    def test_evaluate_precursor(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Mock the create_demand method to return a specific value
        mock_reaction = MagicMock()
        mock_reaction.id = 'demand_metabolite1'
        model.create_demand = MagicMock(return_value=mock_reaction)

        # Mock the slim_optimize method to return a specific value
        model.slim_optimize = MagicMock(return_value=1.0)

        # Create a Metabolite object
        metabolite = Metabolite('metabolite1', compartment='compartment1')

        # Set the reactions attribute to a DictList object
        model.reactions = DictList([mock_reaction])

        # Call the evaluate_precursor method
        flux_value = model.evaluate_precursor(model, metabolite)

        # Check that the evaluate_precursor method returns the correct value
        self.assertEqual(flux_value, 1.0)

    import os
    from cobra import Metabolite

    def test_to_sbml(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Set the seeds and targets attributes to a list of tuples
        model.seeds = [('metabolite1', 'compartment1')]
        model.targets = [('metabolite2', 'compartment2')]

        # Define the file name, save path, and seeds and targets options
        file_name = 'test_model.xml'
        save_path = '/Users/josediogomoura/gap_filling_dl/tests/data/test_model/output'
        seeds = True
        targets = True

        # Call the to_sbml method
        model.to_sbml(file_name, save_path, seeds, targets)

        # Check that the SBML files were created
        self.assertTrue(os.path.exists(os.path.join(save_path, file_name.replace('.xml', '_seeds.xml'))))
        self.assertTrue(os.path.exists(os.path.join(save_path, file_name.replace('.xml', '_targets.xml'))))

        # Clean up the created files
        os.remove(os.path.join(save_path, file_name.replace('.xml', '_seeds.xml')))
        os.remove(os.path.join(save_path, file_name.replace('.xml', '_targets.xml')))

    from unittest.mock import patch

    def test_create_random_knockout_models(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Define the parameters for the create_random_knockout_models method
        num_models = 5
        min_knockouts = 1
        max_knockouts = 3
        generate_files = False
        parent_folder = None
        seeds_targets = False
        universal_model = None

        # Mock the _initialize_universal_model_reactions_ids, _generate_random_knockout_models, _generate_files, and _print_removed_reactions methods
        with patch.object(model, '_initialize_universal_model_reactions_ids', return_value=(None, None, None)):
            with patch.object(model, '_generate_random_knockout_models', return_value=(None, [1, 2, 3, 2, 1], None)):
                with patch.object(model, '_generate_files'):
                    with patch.object(model, '_print_removed_reactions'):
                        # Call the create_random_knockout_models method
                        knockout_numbers = model.create_random_knockout_models(num_models, min_knockouts, max_knockouts,
                                                                               generate_files, parent_folder,
                                                                               seeds_targets, universal_model)

        # Check that the create_random_knockout_models method returns the correct value
        self.assertEqual(knockout_numbers, [1, 2, 3, 2, 1])

    from unittest.mock import patch

    def test_initialize_universal_model_reactions_ids(self):
        # Initialize the Model
        model = Model(self.cobra_model, self.objective_function_id)

        # Define the universal model and save path
        universal_model = 'path_to_your_universal_model'
        save_path = 'path_to_save_directory'

        # Mock the read_sbml_model and write_sbml_model functions to avoid reading and writing actual files
        with patch('cobra.io.read_sbml_model') as mock_read_sbml_model, patch(
                'cobra.io.write_sbml_model') as mock_write_sbml_model:
            # Define the mock universal model
            mock_universal_model = MagicMock()
            mock_universal_model.reactions = [MagicMock(id='reaction1_c'), MagicMock(id='reaction2_c')]
            mock_read_sbml_model.return_value = mock_universal_model

            # Call the _initialize_universal_model_reactions_ids method
            universal_model_reactions_ids, kegg_reactions, universal_model_file = model._initialize_universal_model_reactions_ids(
                universal_model, save_path)

        # Check that the _initialize_universal_model_reactions_ids method returns the correct value
        self.assertEqual(universal_model_reactions_ids, {'reaction1', 'reaction2'})
        self.assertEqual(kegg_reactions, [])
        self.assertEqual(universal_model_file, os.path.join(save_path, 'universal_model.xml'))

    # def test_generate_random_knockout_models(self):
    #     # Create a mock reaction
    #     reaction1 = create_autospec(Reaction, instance=True)
    #     reaction1.id = 'reaction1'
    #     reaction1.reverse_id = 'reaction1_reverse'
    #     reaction2 = create_autospec(Reaction, instance=True)
    #     reaction2.id = 'reaction2'
    #     reaction2.reverse_id = 'reaction2_reverse'
    #     reaction3 = create_autospec(Reaction, instance=True)
    #     reaction3.id = 'reaction3'
    #     reaction3.reverse_id = 'reaction3_reverse'
    #
    #     # Initialize a cobra.Model with mock reactions
    #     cobra_model = cobra.Model('mock_model')
    #     cobra_model.add_reactions([reaction1, reaction2, reaction3])
    #
    #     # Initialize the Model with the cobra.Model
    #     model = Model(cobra_model)
    #
    #     # Define the parameters for the _generate_random_knockout_models method
    #     num_models = 2
    #     min_knockouts = 1
    #     max_knockouts = 2
    #     kegg_reactions = ['reaction1', 'reaction2', 'reaction3']
    #
    #     # Mock the identify_seeds, identify_targets, and remove_reactions methods to return specific values
    #     attempts = [0]
    #     with patch.object(Model, 'identify_seeds',
    #                       side_effect=lambda: ['metabolite1'] if attempts[0] % 2 == 0 else []), \
    #             patch.object(Model, 'identify_targets',
    #                          side_effect=lambda: ['metabolite2'] if attempts[0] % 2 == 0 else []), \
    #             patch.object(Model, 'remove_reactions'):
    #         # Call the _generate_random_knockout_models method
    #         print("Calling _generate_random_knockout_models method...")
    #         knockout_models, knockout_numbers, removed_reactions_dict = model._generate_random_knockout_models(
    #             num_models,
    #             min_knockouts,
    #             max_knockouts,
    #             kegg_reactions
    #         )
    #         print("Finished calling _generate_random_knockout_models method")
    #
    #         # Check that the _generate_random_knockout_models method returns the correct value
    #     print("Checking the results...")
    #     self.assertEqual(len(knockout_models), num_models)
    #     self.assertTrue(all(min_knockouts <= num <= max_knockouts for num in knockout_numbers))
    #     self.assertEqual(len(removed_reactions_dict), num_models)
    #     print("Test passed")

    def test_generate_files(self):
        # Create a real Model instance
        cobra_model = cobra.Model('mock_model')
        model = Model(cobra_model)
        model.id = 'mock_model'
        model.objective_function_id = 'objective_function'

        # Define the parameters for the _generate_files method
        knockout_models = [model]
        parent_folder = '/path/to/parent_folder'
        seeds_targets = True
        universal_model_file = '/path/to/universal_model_file'

        # Mock the os.makedirs, cobra.io.write_sbml_model, Model.identify_seeds, Model.identify_targets, Model.to_sbml, and shutil.copy2 methods
        with patch.object(os, 'makedirs') as mock_makedirs, \
                patch.object(cobra.io, 'write_sbml_model') as mock_write_sbml_model, \
                patch.object(Model, 'identify_seeds') as mock_identify_seeds, \
                patch.object(Model, 'identify_targets') as mock_identify_targets, \
                patch.object(Model, 'to_sbml') as mock_to_sbml, \
                patch.object(shutil, 'copy2') as mock_copy2:
            # Call the _generate_files method
            model._generate_files(knockout_models, parent_folder, seeds_targets, universal_model_file)

        # Assert that the mocked methods were called with the expected arguments
        mock_makedirs.assert_called_once_with(os.path.join(parent_folder, 'model_1'), exist_ok=True)
        mock_write_sbml_model.assert_called_once_with(model, os.path.join(parent_folder, 'model_1', 'mock_model.xml'))
        mock_identify_seeds.assert_called_once()
        mock_identify_targets.assert_called_once()
        mock_to_sbml.assert_called_once_with('mock_model.xml', save_path=os.path.join(parent_folder, 'model_1'),
                                             seeds=True, targets=True)
        mock_copy2.assert_called_once_with(universal_model_file, os.path.join(parent_folder, 'model_1',
                                                                              os.path.basename(universal_model_file)))


    def test_print_removed_reactions(self):
        # Create a real Model instance
        cobra_model = cobra.Model('mock_model')
        model = Model(cobra_model)
        model.id = 'mock_model'
        model.objective_function_id = 'objective_function'

        # Define the removed_reactions_dict variable
        removed_reactions_dict = {0: ['reaction1', 'reaction2'], 1: ['reaction3', 'reaction4']}

        # Mock the print function
        with patch('builtins.print') as mock_print:
            # Call the _print_removed_reactions method
            model._print_removed_reactions(removed_reactions_dict)

        # Assert that the print function was called with the expected arguments
        mock_print.assert_any_call("Reactions removed from each model:")
        mock_print.assert_any_call("Model 1: ['reaction1', 'reaction2']")
        mock_print.assert_any_call("Model 2: ['reaction3', 'reaction4']")

    def test_get_bio_precursors(self):
        # Create a mock Model instance
        mock_model = create_autospec(Model, instance=True)
        mock_model.id = 'mock_model'
        mock_model.objective_function_id = 'objective_function'

        # Create mock Metabolite instances
        metabolite1 = create_autospec(Metabolite, instance=True)
        metabolite2 = create_autospec(Metabolite, instance=True)

        # Create a mock Reaction instance
        reaction = create_autospec(Reaction, instance=True)
        reaction.reactants = [metabolite1, metabolite2]

        # define bio_precursors attribute
        mock_model.bio_precursors = [metabolite1, metabolite2]

        # Mock the get_reactants method to return the mock Reaction instance
        mock_model.get_reactants = lambda \
            objective_function_id: reaction.reactants if objective_function_id == 'objective_function' else []

        # Mock the get_bio_precursors method to return the reactants of the objective function
        mock_model.get_bio_precursors = lambda: mock_model.get_reactants(mock_model.objective_function_id)

        # Call the get_bio_precursors method
        bio_precursors = mock_model.get_bio_precursors()

        # Assert that the get_bio_precursors method returns the correct value
        assert bio_precursors == [metabolite1, metabolite2]
        assert mock_model.bio_precursors == [metabolite1, metabolite2]

    def test_get_pre_precursors(self):
        # Create a mock Model instance
        mock_model = create_autospec(Model, instance=True)
        mock_model.id = 'mock_model'
        mock_model.objective_function_id = 'objective_function'

        # Create mock Metabolite instances
        metabolite1 = create_autospec(Metabolite, instance=True)
        metabolite1.id = 'metabolite1'
        metabolite2 = create_autospec(Metabolite, instance=True)
        metabolite2.id = 'metabolite2'

        # Define the bio_precursors variable and set it as an attribute of the mock Model instance
        bio_precursors = [metabolite1, metabolite2]
        mock_model.bio_precursors = bio_precursors

        # Mock the get_bio_precursors and get_reactants methods to return specific values
        mock_model.get_bio_precursors = lambda: bio_precursors
        mock_model.get_reactants = lambda reaction_id: [metabolite1] if reaction_id != 'objective_function' else []

        # Mock the get_pre_precursors method to return the reactants of the reactions in which the biomass precursors participate
        mock_model.get_pre_precursors = lambda: {'metabolite1': [metabolite1], 'metabolite2': [metabolite1]}

        # Call the get_pre_precursors method
        pre_precursors = mock_model.get_pre_precursors()

        # Assert that the get_pre_precursors method returns the correct value
        assert pre_precursors == {'metabolite1': [metabolite1], 'metabolite2': [metabolite1]}

    def test_get_reaction(self):
        # Create a mock Model instance
        mock_model = create_autospec(Model, instance=True)

        # Define the reaction variable
        reaction = 'reaction1'

        # Create a mock Reaction instance
        mock_reaction = create_autospec(Reaction, instance=True)

        # Create a mock DictList instance for the reactions attribute of the mock Model instance
        mock_reactions = create_autospec(DictList, instance=True)
        mock_reactions.get_by_id = lambda r: mock_reaction if r == reaction else None

        # Set the reactions attribute on the mock Model instance
        mock_model.reactions = mock_reactions

        # Mock the get_reaction method to return the reaction from the model by its ID
        mock_model.get_reaction = lambda r: mock_reactions.get_by_id(r)

        # Call the get_reaction method
        result = mock_model.get_reaction(reaction)

        # Assert that the get_reaction method returns the correct value
        assert result == mock_reaction

    def test_get_products(self):
        model = Model()
        reaction = Reaction('R1')
        metabolite1 = Metabolite('M1')
        metabolite2 = Metabolite('M2')
        reaction.add_metabolites({metabolite1: -1, metabolite2: 1})
        model.add_reactions([reaction])

        expected_products = [metabolite2]
        actual_products = model.get_products('R1')
        self.assertEqual(expected_products, actual_products)

    def test_create_trnas_reactions(self):
        model = Model()
        reaction = Reaction('R1')
        metabolite1 = Metabolite('H2O')
        metabolite2 = Metabolite('e-Protein')
        metabolite3 = Metabolite('M1')
        reaction.add_metabolites({metabolite1: -1, metabolite2: -1, metabolite3: 1})
        model.add_reactions([reaction])
        model.precursors_reactions = {'e-Protein': ['R1']}

        model.create_trnas_reactions()

        # Check if the sink reaction for the metabolite M1 is created
        self.assertIsNotNone(model.get_reaction('Sk_M1'))

    def test_create_demand(self):
        model = Model()
        metabolite = Metabolite('M1')
        model.add_metabolites([metabolite])

        actual_reaction = model.create_demand('M1')
        expected_reaction = model.get_reaction('DM_M1')

        self.assertIsNotNone(actual_reaction)
        self.assertEqual(expected_reaction.id, actual_reaction.id)
        self.assertEqual(expected_reaction.bounds, actual_reaction.bounds)
        self.assertEqual(expected_reaction.metabolites, actual_reaction.metabolites)

    def test_create_sink(self):
        model = Model()
        metabolite = Metabolite('M1')
        model.add_metabolites([metabolite])

        actual_reaction = model.create_sink('M1')
        expected_reaction = model.get_reaction('Sk_M1')

        self.assertIsNotNone(actual_reaction)
        self.assertEqual(expected_reaction.id, actual_reaction.id)
        self.assertEqual(expected_reaction.bounds, actual_reaction.bounds)
        self.assertEqual(expected_reaction.metabolites, actual_reaction.metabolites)

    def test_create_reaction(self):
        model = Model()
        reaction_id = 'R1'

        actual_reaction = model.create_reaction(reaction_id)
        expected_reaction = model.get_reaction('R1')

        self.assertIsNotNone(actual_reaction)
        self.assertEqual(expected_reaction.id, actual_reaction.id)

    def test_get_metabolite(self):
        model = Model()
        metabolite = Metabolite('M1')
        model.add_metabolites([metabolite])

        expected_metabolite = metabolite
        actual_metabolite = model.get_metabolite('M1')
        self.assertEqual(expected_metabolite, actual_metabolite)

    def test_get_reactants(self):
        model = Model()
        reaction = Reaction('R1')
        metabolite1 = Metabolite('M1')
        metabolite2 = Metabolite('M2')
        reaction.add_metabolites({metabolite1: -1, metabolite2: 1})
        model.add_reactions([reaction])

        expected_reactants = [metabolite1]
        actual_reactants = model.get_reactants('R1')
        self.assertEqual(expected_reactants, actual_reactants)



if __name__ == "__main__":
    unittest.main()
