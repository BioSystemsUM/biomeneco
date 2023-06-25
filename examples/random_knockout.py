import os
import time

from src.gap_filling_dl.biomeneco.gapfiller import GapFiller


def print_gapfilling_stats(model_name, results):
    reactions_removed = {
        "model_1": ['R00939_C3_in', 'R03191_C3_in', 'R03066_C3_in', 'R04230_C3_in'],
        "model_2": ['R04440_C3_in', 'R04960_C3_in', 'R01466_C3_in', 'R04954_C3_in', 'R02016_C3_in'],
        "model_3": ['R01231_C3_in', 'R01230_C3_in', 'R01528_C3_in', 'R02749_C3_in'],
        "model_4": ['R01819_C3_in'],
        "model_5": ['R00694_C3_in', 'R00566_C3_in', 'R00228_C3_in', 'R00619_C3_in'],
        "model_6": ['R02739_C3_in', 'R03346_C3_in', 'R05070_C3_in', 'R01863_C3_in', 'R00802_C3_in'],
    }
    print(f"--- Gap-filling statistics for {model_name} ---")
    print(results)
    removed_reactions = reactions_removed[model_name]
    num_removed_reactions = len(removed_reactions)
    num_reconstructable = len(results["Reconstructable targets"])

    reconstructable_reactions = set(removed_reactions) & set(results["Reconstructable targets"])
    num_reconstructable_removed = len(reconstructable_reactions)

    print(f"Number of removed reactions: {num_removed_reactions}")
    print(f"Number of reconstructable targets: {num_reconstructable}")
    print(f"Number of removed reactions that are reconstructable: {num_reconstructable_removed}")
    print(f"Reconstructable reactions that were removed: {reconstructable_reactions}")
    print("----------------------------------------------------")


def performance_meneco_folders(models_folder, objective_function_id, temporary_universal_model):
    for model_name in sorted(os.listdir(models_folder)):
        # list all the files in the directory
        try:

            if not model_name.startswith('model_'):
                continue
            model_folder = os.path.join(models_folder, model_name)
            if os.path.isdir(model_folder):
                print('model_folder: ', model_folder)

                print(f"Running GapFiller for {model_name}")
                gf1 = GapFiller.from_folder(model_folder, objective_function_id='Biomass_assembly_C3_in',
                                            temporary_universal_model=temporary_universal_model)  # mudar aqui para ser mais rapido

                start_time = time.time()
                results = gf1.run(enumeration=True, json_output=True)
                end_time = time.time()
                elapsed_time = end_time - start_time

                print(f"Time taken for {model_name}: {elapsed_time}")
                print(f"Results for {model_name}:")
                print(results)
                print_gapfilling_stats(model_name, results)
            else:
                print(f"{model_name} is not a directory. Skipping...")

        except Exception as e:
            print(f"Error for {model_name}", str(e))
            continue


if __name__ == '__main__':
    models_folder = '/Users/josediogomoura/gap_filling_dl/tests/data/original_model'
    objective_function_id = 'Biomass_assembly_C3_in'
    temporary_universal_model = True
    performance_meneco_folders(models_folder, objective_function_id, temporary_universal_model)
