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

    def test_performance_gapfill(self):
        initial_time = time.time()
        reactions_to_add = gapfill(self.my_model, self.universal_model, demand_reactions=True, exchange_reactions=True)
        final_time = time.time()
        elapsed_time = final_time - initial_time
        print(f"Elapsed time: {elapsed_time} seconds")
        print(reactions_to_add)
