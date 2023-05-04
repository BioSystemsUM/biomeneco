import fnmatch
import os
import time
from unittest import TestCase
from cobra.flux_analysis import gapfill
from cobra.io import read_sbml_model
from tests import TEST_DIR


class TestCobraGapfill(TestCase):
    def setUp(self) -> None:
        self.my_model = read_sbml_model(os.path.join(TEST_DIR, "data/6gaps/6gaps_model.xml"))
        self.my_model.objective = 'Biomass_C3_in'
        self.universal_model = read_sbml_model(os.path.join(TEST_DIR, "data/6gaps/universal_model_kegg.xml"))

        # folder with the models to test
        self.models_folder = os.path.join(TEST_DIR, "data/original_model")
        self.removed_reactions = {
            "model_1": ['R00939_C3_in', 'R03191_C3_in', 'R03066_C3_in', 'R04230_C3_in'],
            "model_2": ['R04440_C3_in', 'R04960_C3_in', 'R01466_C3_in', 'R04954_C3_in', 'R02016_C3_in'],
            "model_3": ['R01231_C3_in', 'R01230_C3_in', 'R01528_C3_in', 'R02749_C3_in'],
            "model_4": ['R01819_C3_in'],
            "model_5": ['R00694_C3_in', 'R00566_C3_in', 'R00228_C3_in', 'R00619_C3_in'],
            "model_6": ['R02739_C3_in', 'R03346_C3_in', 'R05070_C3_in', 'R01863_C3_in', 'R00802_C3_in'],
        }

    def test_performance_gapfill(self):
        initial_time = time.time()
        reactions_to_add = gapfill(self.my_model, self.universal_model, demand_reactions=True, exchange_reactions=True)
        final_time = time.time()
        elapsed_time = final_time - initial_time
        print(f"Elapsed time: {elapsed_time} seconds")
        print(reactions_to_add)

    def test_performance_gapfill_folders(self):
        for model_name in sorted(os.listdir(self.models_folder)):
            model_folder = os.path.join(self.models_folder, model_name)
            if os.path.isdir(model_folder):
                print(f"Running Cobra gapfill for {model_name}")

                model_file = None
                for file in os.listdir(model_folder):
                    if fnmatch.fnmatch(file, "model_spneumoniaeR6_knockout_*.xml"):
                        model_file = file
                        break

                if model_file:
                    my_model = read_sbml_model(os.path.join(model_folder, model_file))
                    my_model.objective = 'Biomass_assembly_C3_in'
                    universal_model = read_sbml_model(os.path.join(model_folder, "universal_model.xml"))

                    initial_time = time.time()
                    reactions_to_add = gapfill(my_model, universal_model, demand_reactions=True, exchange_reactions=True)
                    final_time = time.time()
                    elapsed_time = final_time - initial_time

                    print(f"Elapsed time for {model_name}: {elapsed_time} seconds")
                    print(f"Reactions to add for {model_name}:")
                    print(reactions_to_add)

                    self.print_gapfilling_stats(model_name, reactions_to_add)
                else:
                    print(f"No model file found in {model_name} folder")

    def print_gapfilling_stats(self, model_name, reactions_to_add):
        print(f"--- Gap-filling statistics for {model_name} ---")
        removed_reactions = self.removed_reactions[model_name]
        num_removed_reactions = len(removed_reactions)

        reconstructable_reactions = set(removed_reactions) & set(reactions_to_add[0])
        num_reconstructable_removed = len(reconstructable_reactions)

        print(f"Number of removed reactions: {num_removed_reactions}")
        print(f"Number of removed reactions that are reconstructable: {num_reconstructable_removed}")
        print(f"Reconstructable removed reactions:")
        print(reconstructable_reactions)
        print("-------------------------------------------------")

