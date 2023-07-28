import sys
import time
sys.path.append('../..')
sys.path.append('../../src')
sys.path.append('../../examples/BioISO/src')
import cobra
from gap_filling_dl.biomeneco.model import Model
from src.gap_filling_dl.biomeneco.gapfiller import GapFiller
import logging
import os
import warnings
warnings.filterwarnings("ignore")
logging.getLogger("cobra").setLevel(logging.CRITICAL)
logging.getLogger("meneco").setLevel(logging.CRITICAL)


if __name__ == '__main__':
    path_model = '../data/draft_synechocystis/input/model.xml'
    if os.path.exists('../data/draft_synechocystis/universal_model_compartmentalized.xml'):
        os.remove('../data/draft_synechocystis/universal_model_compartmentalized.xml')
    if os.path.exists('../data/draft_synechocystis/temporary_universal_model.xml'):
        os.remove('../data/draft_synechocystis/temporary_universal_model.xml')
    cobra_model = cobra.io.read_sbml_model(path_model)
    my_model = Model(cobra_model, "e_Biomass__cytop")
    # my_model.create_sink("C00002__cytop")
    # my_model.create_sink("C00004__cytop")
    my_model.identify_seeds()
    my_model.identify_targets()
    my_model.create_trnas_reactions()
    my_model.to_sbml('draft_synechocystis', '../data/draft_synechocystis', seeds=True, targets=True)
    initial_time = time.time()
    gf = GapFiller.from_folder('../data/draft_synechocystis/', temporary_universal_model=True, objective_function_id="e_Biomass__cytop", compartments=["cytop", "extr"])
    print("Running gapfilling...")
    res = gf.run(False, write_to_file=True, removed_reactions=None, objective_function_id="e_Biomass__cytop")
    print(res)
    print("Time elapsed: ", time.time() - initial_time)
