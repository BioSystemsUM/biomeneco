import json
import platform
import shutil
from logging.handlers import TimedRotatingFileHandler
from os.path import join
import cobra
from biomeneco.model import Model
from biomeneco.gapfiller import GapFiller
import sys
import os
import logging
import warnings


warnings.filterwarnings("ignore")
logging.getLogger("cobra").setLevel(logging.CRITICAL)
logging.getLogger("meneco").setLevel(logging.CRITICAL)


def read_model_and_run_gap_filling(processing_path, results_path):
    # Read the model with cobra
    model_path = join(processing_path, "model.xml")
    model = cobra.io.read_sbml_model(model_path)
    submission_parameters_path = join(processing_path, "SubmissionParameters.json")
    parameters = json.load(open(submission_parameters_path))

    objective_function_id = parameters["objective_function_id"]

    # Create seeds + targets from the draft model
    my_model = Model(model, objective_function_id=objective_function_id)
    my_model.create_trnas_reactions()  # Create tRNAs reactions for protein synthesis

    for sink in parameters["sinks"]:
        my_model.create_sink(sink)

    my_model.identify_seeds()
    my_model.identify_targets()

    # Create the sbml files
    my_model.to_sbml('model', processing_path, targets=True, seeds=True)

    # Create the GapFiller object
    gf = GapFiller.from_folder(processing_path, results_path, temporary_universal_model=True,
                               objective_function_id=objective_function_id, compartments=parameters["compartments"])

    # Run the gap filling (optimize=False for the fastest method)
    gf.run(optimize=False, write_to_file=True)


if __name__ == '__main__':
    # Call the function to run the commands
    #

    if platform.system() == 'Linux':
        processingPath = sys.argv[1]
        # resultsPath = sys.argv[2]
    elif platform.system() == 'Windows':
        processingPath = r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\lactis\input"
        resultsPath = r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\lactis\output"
    else:
        print('Running on another operating system')

    print("processingPath: ", processingPath)
    print("resultsPath: ", resultsPath)
    logPath = resultsPath + "/trace_errors.log"

    if not os.path.exists(resultsPath):
        os.makedirs(resultsPath)
    else:
        shutil.rmtree(resultsPath)
        os.makedirs(resultsPath)

    # format the log entries
    formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
    handler = TimedRotatingFileHandler(logPath, when='midnight', backupCount=20)
    handler.setFormatter(formatter)
    logger = logging.getLogger(__name__)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    read_model_and_run_gap_filling(processingPath, resultsPath)
