from typing import List, Tuple
import cobra import Model, Metabolite
from cobra.io import write_sbml_model
import os


def write_metabolites_to_sbml(file_name: str, save_path: str, metabolites: List[Tuple[str, str]]):
    """
    Write a list of metabolites to an SBML file.
    Parameters
    ----------
    model_path: str
        The path to the model to write the metabolites to.
    file_name: str
        The name of the SBML file to write to.
    save_path: str
        The path to save the SBML file to.
    metabolites: Tuple[str, str]
        A list of metabolites to write to the SBML file.
    """

    model = cobra.Model()
    metabolite_cobra_metabolite_objects = []
    for metabolite_id, compartment in metabolites:
        metabolite = Metabolite(id=metabolite_id, compartment=compartment)
        metabolite_cobra_metabolite_objects.append(metabolite)

    model.add_metabolites(metabolite_cobra_metabolite_objects)
    write_sbml_model(model, os.path.join(save_path, file_name))
