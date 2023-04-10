def get_metabolite_pathway_map(model):
    """
    Get a map of metabolite IDs to their associated pathways.
    Returns a dictionary where the keys are metabolite IDs and the values are lists of pathways.
    -------
    """
    if not model.reaction_pathway_map:
        model.reaction_pathway_map = get_reaction_pathway_map(model)
    metabolite_pathway_map = {}
    for metabolite in model.metabolites:
        for reaction in metabolite.reactions:
            if reaction.id in model.reaction_pathway_map.keys():
                for pathway in model.reaction_pathway_map[reaction.id]:
                    if pathway not in metabolite_pathway_map.get(metabolite.id, []):
                        metabolite_pathway_map[metabolite.id] = metabolite_pathway_map.get(metabolite.id, []) + [pathway]
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

