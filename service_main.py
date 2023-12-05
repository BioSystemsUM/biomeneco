import cobra
from src.gap_filling_dl.biomeneco.model import Model
from src.gap_filling_dl.biomeneco.gapfiller import GapFiller
import sys
# print("Python Executable:", sys.executable)
# print("Python Path:", sys.path)
import os
import logging
import warnings

warnings.filterwarnings("ignore")
logging.getLogger("cobra").setLevel(logging.CRITICAL)
logging.getLogger("meneco").setLevel(logging.CRITICAL)


def read_model_and_run_gap_filling():
    # Read the model with cobra
    if os.path.exists('../data/draft_synechocystis/universal_model_compartmentalized.xml'):
        os.remove('../data/draft_synechocystis/universal_model_compartmentalized.xml')
    if os.path.exists('../data/draft_synechocystis/temporary_universal_model.xml'):
        os.remove('../data/draft_synechocystis/temporary_universal_model.xml')

    model_path = 'tests/data/draft_synechocystis/input/model.xml'
    model = cobra.io.read_sbml_model(model_path)

    # Define the objective function id
    objective_function_id = "e_Biomass__cytop"

    # Create seeds + targets from the draft model
    my_model = Model(model, objective_function_id=objective_function_id)
    my_model.create_trnas_reactions()  # Create tRNAs reactions for protein synthesis

    my_model.identify_seeds()
    my_model.identify_targets()

    # Create the sbml files
    sbml_output_path = 'tests/data/draft_synechocystis/input'
    my_model.to_sbml('model_synechocystis', sbml_output_path, targets=True, seeds=True)

    # Create the GapFiller object
    gf = GapFiller.from_folder(sbml_output_path, 'tests/data/draft_synechocystis/output', temporary_universal_model=True,
                               objective_function_id=objective_function_id, clone=True)

    # Run the gap filling (optimize=False for the fastest method)
    gf.run(optimize=False, write_to_file=True)


if __name__ == '__main__':
    # Call the function to run the commands
    read_model_and_run_gap_filling()