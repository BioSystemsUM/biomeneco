from datetime import time
from unittest import TestCase
from cobra.flux_analysis import gapfill
from cobra.io import read_sbml_model


class TestCobraGapfill(TestCase):
    def setUp(self) -> None:
        self.my_model = read_sbml_model('../data/new_toy_test.xml')
        self.my_model.objective = 'e_Biomass__in'
        self.universal_model = read_sbml_model('../data/test_model/universal_model_kegg.xml')

    def test_performance_gapfill(self):
        reactions_to_add = gapfill(self.my_model, self.universal_model, demand_reactions=True, exchange_reactions=True)
        print(reactions_to_add)
