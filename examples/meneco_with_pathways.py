import os

from cobra.flux_analysis import gapfill
from cobra.io import read_sbml_model, write_sbml_model
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append('/home/examples/BioISO/src')
from gap_filling_dl.kegg_api import get_related_pathways
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append('/home/examples/BioISO/src')

from gap_filling_dl.biomeneco.utils import get_reaction_pathway_map
from gap_filling_dl.biomeneco.model import Model
from gap_filling_dl.biomeneco.gapfiller import GapFiller


def pre_processing():
    model = read_sbml_model('data/6gaps/6gaps_model.xml')
    # model.objective = 'Biomass_assembly_C3_in'
    # model.add_boundary(model.metabolites.get_by_id('C00004_in'), type='sink')
    # print(model.optimize().objective_value)
    my_model = Model(model=model, objective_function_id='Biomass_assembly_C3_in')
    # my_model.remove_reactions(["R00200_in", "R00703_in"])
    my_model.identify_targets()
    my_model.identify_seeds()
    my_model.to_sbml('toy', 'data/6gaps', seeds=True, targets=True)
    pathways_to_ignore = []  # propably we will have to ignore pathways like 'Metabolic pathways', because they are huge, and include them will make this useless
    pathways_to_keep = []
    for target in my_model.targets:
        pathways_to_keep += [pathway for pathway in my_model.metabolite_pathway_map[target[0]]]
    pathways_to_keep = list(set(pathways_to_keep))
    related_pathways = set()
    for pathway in pathways_to_keep:
        related_pathways.update(get_related_pathways(pathway))
    pathways_to_keep += list(related_pathways)
    print(pathways_to_keep)
    universal_model = read_sbml_model(r"data/universal_model_kegg.xml")
    universal_model = get_reaction_pathway_map(universal_model)
    print(len(universal_model.reactions))
    to_keep = set()
    for pathway in universal_model.groups:
        if pathway.name in pathways_to_keep:
            to_keep.update(reaction for reaction in pathway.members)
    # to_keep.update(reaction for reaction in universal_model.reactions if not universal_model.reaction_pathway_map[reaction.id]) # always keep reactions that are not in any pathway
    to_remove = set(universal_model.reactions) - to_keep
    universal_model.remove_reactions(list(to_remove), remove_orphans=True)
    print(len(universal_model.reactions))
    write_sbml_model(universal_model, r"data/6gaps/universal_model_kegg_temp.xml")


def run_meneco():
    gapfiller = GapFiller.from_folder('data/6gaps')
    result = gapfiller.run(enumeration=True, json_output=True)
    print(result)


def run_cobra_gapfill():
    my_model = read_sbml_model('data/6gaps/6gaps_model.xml')
    my_model.objective = 'Biomass_assembly_C3_in'
    universal_model = read_sbml_model(r"data/6gaps/universal_model_kegg_temp.xml")
    reactions_to_add = gapfill(my_model, universal_model, demand_reactions=True, exchange_reactions=True, iterations=10)
    print(reactions_to_add)



import os
import unittest
import time
from src.gap_filling_dl.biomeneco.gapfiller import GapFiller
from tests import TEST_DIR


def print_or_write(text, file=None):
    if file:
        file.write(text + "\n")
    else:
        print(text)


class TestGapFiller():
    def setUp(self):
        # self.gf = GapFiller.from_folder(os.path.join(TEST_DIR, "data/6gaps"))

        self.reactions_removed = {
            'model_1': ['R00161_C3_in', 'R02023_C3_in', 'R01021_C3_in'],
            'model_2': ['R01416_C3_in'],
            'model_3': ['R00200_C3_in', 'R00462_C3_in'],
            'model_4': ['R00248_C3_in', 'R02689_C3_in'],
            'model_5': ['R01978_C3_in'],
            'model_6': ['R03236_C3_in', 'R00549_C3_in', 'R04364_C3_in', 'R01134_C3_in', 'R01665_C3_in', 'R03066_C3_in']
        }

        self.models_folder = 'data/original_model'

    def test_performance_meneco(self):
        start_time = time.time()
        results = self.gf.run(enumeration=False, json_output=True)
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

    def test_performance_meneco_folders(self):
        self.results = {}

        for model_name in sorted(os.listdir(self.models_folder)):
            # list all the files in the directory
            # try:

                if not model_name.startswith('model_'):
                    continue
                model_folder = os.path.join(self.models_folder, model_name)
                if os.path.isdir(model_folder):
                    print('model_folder: ', model_folder)
                    print(f"Running GapFiller for {model_name}")
                    gf1 = GapFiller.from_folder(model_folder, objective_function_id='Biomass_assembly_C3_in',
                                                temporary_universal_model=True)  # mudar aqui para ser mais rapido

                    start_time = time.time()
                    results = gf1.run(enumeration=False, json_output=True)
                    self.results[model_name] = results
                    end_time = time.time()
                    elapsed_time = end_time - start_time

                    print(f"Time taken for {model_name}: {elapsed_time}")
                    print(f"Results for {model_name}:")
                    print(results)
                    self.print_gapfilling_stats(model_name, results)
                else:
                    print(f"{model_name} is not a directory. Skipping...")

            # except Exception as e:
            #     print(f"Error for {model_name}: {type(e).__name__}: {str(e)}")
            #     continue

        # Save results to a file after all models have been processed
        self.save_results_to_file()

    def print_gapfilling_stats(self, model_name, results, file=None):
        # print_or_write("Reactions removed dictionary keys: " + ", ".join(self.reactions_removed.keys()), file)

        # if model_name not in self.reactions_removed:
        #     print_or_write(f"{model_name} not found in reactions_removed dictionary.", file)
        #     return

        # print_or_write(f"--- Gap-filling statistics for {model_name} ---", file)
        removed_reactions = self.reactions_removed[model_name]
        num_removed_reactions = len(removed_reactions)
        num_reconstructable = len(results["Reconstructable targets"])

        reconstructable_reactions = set(removed_reactions) & set(results["Reconstructable targets"])
        num_reconstructable_removed = len(reconstructable_reactions)

        # print_or_write(f"Number of removed reactions: {num_removed_reactions}", file)
        # print_or_write(f"Number of reconstructable targets: {num_reconstructable}", file)
        # print_or_write(f"Number of removed reactions that are reconstructable: {num_reconstructable_removed}", file)
        # print_or_write(f"Reconstructable reactions that were removed: {reconstructable_reactions}", file)
        # print_or_write("----------------------------------------------------", file)

    def save_results_to_file(self):
        output_folder = 'meneco_outputs'
        os.makedirs(output_folder, exist_ok=True)
        output_file = os.path.join(output_folder, "all_models_results.txt")

        with open(output_file, 'w') as f:
            for model_name, results in self.results.items():
                f.write(f"--- Meneco output for {model_name} ---\n")
                f.write(f"Results for {model_name}:\n")
                for key, value in results.items():
                    f.write(f"{key}: {value}\n")

                self.print_gapfilling_stats(model_name, results, file=f)


if __name__ == '__main__':
    # pre_processing()
    # run_meneco()
    # run_cobra_gapfill()
    test = TestGapFiller()
    test.setUp()
    test.test_performance_meneco_folders()

