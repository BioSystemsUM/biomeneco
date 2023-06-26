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


class TestGapFiller(unittest.TestCase):
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

        self.models_folder = '../data/original_model'

    def test_performance_meneco(self):
        start_time = time.time()
        results = self.gf.run(enumeration=True, json_output=True)
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
            try:

                if not model_name.startswith('model_'):
                    continue
                model_folder = os.path.join(self.models_folder, model_name)
                if os.path.isdir(model_folder):
                    print('model_folder: ', model_folder)

                    print(f"Running GapFiller for {model_name}")
                    gf1 = GapFiller.from_folder(model_folder, objective_function_id='Biomass_assembly_C3_in',
                                                temporary_universal_model=True)  # mudar aqui para ser mais rapido

                    start_time = time.time()
                    results = gf1.run(enumeration=True, json_output=True)
                    self.results[model_name] = results
                    end_time = time.time()
                    elapsed_time = end_time - start_time

                    print(f"Time taken for {model_name}: {elapsed_time}")
                    print(f"Results for {model_name}:")
                    print(results)
                    self.print_gapfilling_stats(model_name, results)
                else:
                    print(f"{model_name} is not a directory. Skipping...")

            except Exception as e:
                print(f"Error for {model_name}: {type(e).__name__}: {str(e)}")

                continue

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
    unittest.main()
