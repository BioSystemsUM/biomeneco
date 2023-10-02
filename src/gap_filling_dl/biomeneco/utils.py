import cobra
from cobra.flux_analysis import flux_variability_analysis
from bioservices import KEGG
from typing import Dict

def edit_metabolite_id(metabolite_id):
    # if metabolite has this structure: 'C02205__C_in', it will be replaced by 'C02205__in', we need to remove the 'C_'
    if 'C_' in metabolite_id:
        metabolite_id = metabolite_id.replace('C_', '')
    return metabolite_id



def edit_metabolite_id(metabolite_id):
    # if metabolite has this structure: 'C02205__C_in', it will be replaced by 'C02205__in', we need to remove the 'C_'
    if 'C_' in metabolite_id:
        metabolite_id = metabolite_id.replace('C_', '')
    return metabolite_id


def get_metabolite_pathway_map(model):
    """
    Get a map of metabolite IDs to their associated pathways.
    Returns a dictionary where the keys are metabolite IDs and the values are lists of pathways.
    -------
    """
    # edit the metabolite id
    for metabolite in model.metabolites:
        metabolite.id = edit_metabolite_id(metabolite.id)

    if not model.reaction_pathway_map:
        model.reaction_pathway_map = get_reaction_pathway_map(model)
    metabolite_pathway_map = {}
    for metabolite in model.metabolites:
        for reaction in metabolite.reactions:
            if reaction.id in model.reaction_pathway_map.keys():
                for pathway in model.reaction_pathway_map[reaction.id]:
                    if pathway not in metabolite_pathway_map.get(metabolite.id, []):
                        metabolite_pathway_map[metabolite.id] = metabolite_pathway_map.get(metabolite.id, []) + [
                            pathway]

    # print if the metabolite is not in any pathway
    for metabolite in model.metabolites:
        if metabolite.id not in metabolite_pathway_map.keys():
            print(f"Metabolite {metabolite.id} is not in any pathway")

    return metabolite_pathway_map


def get_reaction_pathway_map(model):
    """
    Get a map of reaction IDs to their associated pathways.
    Returns a dictionary where the keys are reaction IDs and the values are lists of pathways.
    -------
    """
    pathway_map = {}
    for group in model.groups:
        for member in group.members:
            member_id = member.id
            if member_id in pathway_map:
                pathway_map[member_id].append(group.name)
            else:
                pathway_map[member_id] = [group.name]
    for reaction in model.reactions:
        if reaction.id not in pathway_map.keys():
            pathway_map[reaction.id] = []
    for key, value in pathway_map.items():
        pathway_map[key] = list(set(value))
    model.reaction_pathway_map = pathway_map
    return model


def get_kegg_reaction_ids():
    kegg = KEGG()
    kegg_reaction_ids = kegg.parse(kegg.list("reaction"))  # Get a list of KEGG reaction identifiers
    return kegg_reaction_ids


def map_reaction_names_to_kegg_ids(model, chunk_size=100) -> Dict[str, str]:
    kegg = KEGG()
    kegg_id_mapping = {}

    reaction_ids = [reaction.id for reaction in model.reactions]  # Get a list of reaction IDs

    # Split reaction_ids into chunks
    reaction_chunks = [reaction_ids[i:i + chunk_size] for i in range(0, len(reaction_ids), chunk_size)]

    # não está a funcionar
    for reaction_chunk in reaction_chunks:
        for reaction_id in reaction_chunk:
            # Perform a search for the reaction name
            search_result = kegg.find("reaction", reaction_id)
            print(f"Query: {reaction_id}")
            print(f"Search result: {search_result}")
            search_result = search_result.split("\n")

            # Map the found KEGG IDs to reaction names
            for result in search_result:
                if result:
                    kegg_id, reaction_name = result.split("\t")
                    kegg_id_mapping[reaction_name] = kegg_id

    return kegg_id_mapping


def get_compartments_in_model(model: cobra.Model) -> list:
    """
    Get a list of compartments in a model.
    Parameters
    ----------
    model

    Returns
    -------

    """

    return [e.lstrip("C_") for e in model.compartments.keys()]


def identify_dead_end_metabolites(model: cobra.Model) -> list:
    """
    Identify dead-end metabolites in a model.
    Parameters
    ----------
    model

    Returns
    -------

    """

    dead_ends = []

    for metabolite in model.metabolites:

        reactions = list(metabolite.reactions)  # get the reactions in which the metabolite participates

        if len(reactions) == 1:  # filtrar metabolitos que participam apenas numa reação

            fva_result = flux_variability_analysis(model, reaction_list=reactions)

            if all((fva_result.maximum == 0) & (fva_result.minimum == 0)):
                dead_ends.append(metabolite.id)

    print(f"Number of dead-end metabolites found: {len(dead_ends)}")

    return dead_ends
