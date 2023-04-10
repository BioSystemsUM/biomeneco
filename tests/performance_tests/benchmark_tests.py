
import os
import time

from gap_filling_dl.my_classes import Model, GapFiller, write_metabolites_to_sbml

model_path = '/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/models/iDS372.xml'
universal_model_path = '/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/universal_model_kegg.xml'

model = Model(model_path)

seeds_path = '/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/seeds.xml'
targets_path = '/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/targets.xml'

gap_filler = GapFiller(model_path, universal_model_path)
gap_filler.seeds_path = seeds_path
gap_filler.targets_path = targets_path

start = time.time()
result = gap_filler.run(enumeration=True, json_output=True)
end = time.time()

print("Time taken: ", end - start)

gap_filler.evaluate_results(result)









