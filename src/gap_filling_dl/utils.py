import json
import os

from cobra.io import write_sbml_model, read_sbml_model

from gap_filling_dl.kegg_api import get_related_pathways, create_related_pathways_map
from gap_filling_dl.biomeneco.model import Model
# from k_means_constrained import KMeansConstrained


def build_temporary_universal_model(gap_filler, folder_path, related_pathways):
    """
    Build a temporary universal model from the universal model and the gap-filling results.
    Parameters
    ----------
    related_pathways
    gap_filler
    folder_path

    Returns
    -------

    """
    pathways_to_ignore = {'Metabolic pathways', 'Biosynthesis of secondary metabolites', 'Microbial metabolism in diverse environments',
                          'Biosynthesis of cofactors', 'Carbon metabolism', 'Fatty acid metabolism'}
    pathways_to_keep = []
    metabolite_pathways_map = {}
    read_universal_model = read_sbml_model(gap_filler.universal_model_path)
    universal_model = Model(model=read_universal_model)
    print('Number of reactions in universal model:', len(universal_model.reactions))
    targets = read_sbml_model(gap_filler.targets_path)
    cofactors = json.load(open(os.path.join(folder_path, '../cofactors.json'), 'r'))
    for target in targets.metabolites:
        if target.id in universal_model.metabolite_pathway_map.keys():
            if '__'.join(target.id.split("__")[:-1]) in cofactors.keys():
                pathways_to_keep += [pathway for pathway in universal_model.metabolite_pathway_map[target.id] if pathway in cofactors['__'.join(target.id.split("__")[:-1])]]
                metabolite_pathways_map[target.id] = set(pathway for pathway in universal_model.metabolite_pathway_map[target.id] if pathway in cofactors['__'.join(target.id.split("__")[:-1])])
            else:
                pathways_to_keep += [pathway for pathway in universal_model.metabolite_pathway_map[target.id]]
                metabolite_pathways_map[target.id] = set(pathway for pathway in universal_model.metabolite_pathway_map[target.id])
    pathways_to_keep = set(pathways_to_keep) - pathways_to_ignore
    metabolite_pathways_map = {metabolite: pathways - pathways_to_ignore for metabolite, pathways in metabolite_pathways_map.items() if pathways not in pathways_to_ignore}
    if related_pathways:
        related_pathways = set()
        if os.path.exists(os.path.join(folder_path, '../related_pathways_map.json')):
            with open(os.path.join(folder_path, '../related_pathways_map.json'), 'r') as related_pathways_map_file:
                related_pathways_map = json.load(related_pathways_map_file)
        else:
            related_pathways_map = create_related_pathways_map(universal_model, folder_path)

        for pathway in pathways_to_keep:
            related_pathways.update(get_related_pathways(pathway, related_pathways_map))
        pathways_to_keep = pathways_to_keep.union(related_pathways) - pathways_to_ignore

    for metabolite, pathways in metabolite_pathways_map.items():
        copy = pathways.copy()
        for pathway in copy:
            metabolite_pathways_map[metabolite].update(get_related_pathways(pathway, related_pathways_map))
    print('Pathways to keep are:', pathways_to_keep)
    to_keep = set()
    number_of_reactions_by_pathway = {}
    for pathway in universal_model.groups:
        if pathway.name in pathways_to_keep:
            number_of_reactions_by_pathway[pathway.name] = len(pathway.members)
    for pathway in universal_model.groups:
        if pathway.name in pathways_to_keep:
            to_keep.update(reaction for reaction in pathway.members)
    to_remove = set(universal_model.reactions) - to_keep
    # get reactions without pathway
    # for reaction_id in universal_model.reaction_pathway_map.keys():
    #     if len(universal_model.reaction_pathway_map[reaction_id]) == 0 and universal_model.reactions.get_by_id(reaction_id) in to_remove:
    #         to_remove.remove(universal_model.reactions.get_by_id(reaction_id))
    universal_model.remove_reactions(list(to_remove), remove_orphans=True)
    print('Number of reactions in temporary universal model:', len(universal_model.reactions))
    write_sbml_model(universal_model, os.path.join(folder_path, 'temporary_universal_model.xml'))

