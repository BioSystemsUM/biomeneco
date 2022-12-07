from typing import Tuple, List

import cobra
from cobra import Model, Metabolite
from bioiso import BioISO
from bioiso import load, set_solver


def write_metabolites_to_sbml(file_path: str, metabolites: List[Tuple[str, str]]):
    """
    Write a list of metabolites to an SBML file.

    Parameters
    ----------
    file_path: str
        The path to the SBML file to write to.
    metabolites: Tuple[str, str]
        A list of metabolites to write to the SBML file.
    """
    metabolites_sbml = Model(file_path)

    metabolite_cobra_metabolite_objects = []
    for metabolite_id, compartment in metabolites:
        metabolite = Metabolite(id=metabolite_id, compartment=compartment)
        metabolite_cobra_metabolite_objects.append(metabolite)

    metabolites_sbml.add_metabolites(metabolite_cobra_metabolite_objects)
    cobra.io.write_sbml_model(metabolites_sbml, file_path)


def identify_seeds(model_path: str) -> List[Tuple[str, str]]:
    """
    Identify the seeds (external metabolites that enter the model) of a model.

    Parameters
    ----------
    model_path: str
        The path to the model to identify the seeds of.
    """

    model = load(model_path)
    total_seeds = []
    total_seeds_ids = []
    for reaction in model.reactions:
        found_out = False
        found_inside = False
        seeds = []
        seeds_ids = []
        reactants = reaction.reactants
        products = reaction.products
        for reactant in reactants:
            if reactant.compartment == "C_00001":
                found_out = True
                if reaction.lower_bound < 0 and reactant.id not in total_seeds_ids:
                    seeds_ids.append(reactant.id)
                    seeds.append((reactant.id, reactant.compartment))

        for product in products:
            if product.compartment == "C_00002":
                found_inside = True
            elif product.compartment == "C_00001":
                if reaction.upper_bound > 0 and product.id not in total_seeds_ids \
                        and product.id not in seeds_ids:
                    seeds_ids.append(product.id)
                    seeds.append((product.id, product.compartment))

        if found_out and found_inside:
            total_seeds.extend(seeds)
            total_seeds_ids.extend(seeds_ids)

    return total_seeds


def identify_targets(model_path: str, objective_function: str, objective: str, solver: str) -> List[Tuple[str, str]]:
    """
    Optimised way of identifying targets for a model using BioISO.

    Parameters
    ----------
    model_path: str
        The path to the model to identify the targets of.
    objective_function:
        The objective function to use.
    objective:
        The objective to use.
    solver: str
        The solver to use.

    Returns
    -------
    List[Tuple[str, str]]
        A list of targets where the first element is the metabolite target identifier and the second element is the
        compartment identifier.
    """
    model = load(model_path)

    set_solver(model, solver)

    bio = BioISO(objective_function, model, objective)
    bio.run(2, False)

    results = bio.get_tree()
    biomass_components = results["M_root_M_root_M_root_product"]["next"]

    targets = []
    for biomass_component in biomass_components:

        biomass_component_role = biomass_components[biomass_component].get("role")
        if biomass_component_role == "Reactant":
            specific_biomass_components = biomass_components[biomass_component].get("next")
            for specific_biomass_component in specific_biomass_components:

                specific_biomass_component_report = specific_biomass_components[specific_biomass_component]
                analysis = specific_biomass_component_report.get("analysis")
                role = specific_biomass_component_report.get("role")
                if role == "Reactant" and not analysis:
                    target_id = specific_biomass_component_report.get("identifier")
                    compartment = specific_biomass_component_report.get("compartment")

                    if target_id not in targets:
                        targets.append((target_id, compartment))

    return targets

    # print(results)

# identify_targets("../tests")
