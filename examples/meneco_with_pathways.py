import os

from cobra.flux_analysis import gapfill
from cobra.io import read_sbml_model, write_sbml_model

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


if __name__ == '__main__':
    pre_processing()
    run_meneco()
    # run_cobra_gapfill()
