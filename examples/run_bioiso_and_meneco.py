import os

import time
from meneco import run_meneco

from gap_filling_dl.biomeneco.generate_targets import identify_targets, write_metabolites_to_sbml, identify_seeds

EXAMPLES_DIR = os.path.dirname(os.path.abspath(__file__))

objective_function = "Biomass_assembly_C3_cytop"
objective = "maximize"
solver = 'cplex'

model_path = os.path.join(EXAMPLES_DIR, "data", "models", "6gaps_model.xml")
targets = identify_targets(model_path, objective_function,
                           objective, solver)

write_metabolites_to_sbml(os.path.join(EXAMPLES_DIR, "data", "temp_targets.xml"), targets)

seeds = identify_seeds(model_path)

write_metabolites_to_sbml(os.path.join(EXAMPLES_DIR, "data", "temp_seeds.xml"), seeds)
start = time.time()
temp_seeds_path = os.path.join(EXAMPLES_DIR, "data", "temp_seeds.xml")
temp_targets_path = os.path.join(EXAMPLES_DIR, "data", "temp_targets.xml")
temp_universal_model = os.path.join(EXAMPLES_DIR, "data", "kegg_universal_model.xml")

result = run_meneco(draftnet=model_path,
                    seeds=temp_seeds_path,
                    targets=temp_targets_path,
                    repairnet=temp_universal_model,
                    enumeration=True,
                    json_output=True)
end = time.time()
