import os
import random
import shutil
import tempfile
from copy import deepcopy
from typing import Tuple, List, Iterable
import cobra
from bioiso import BioISO
from bioiso import set_solver
from cobra import Reaction

from gap_filling_dl.biomeneco.utils import get_metabolite_pathway_map, get_reaction_pathway_map
from src.gap_filling_dl.write_sbml.metabolites_to_sbml import write_metabolites_to_sbml


class Model(cobra.Model):
    def __init__(self, model: cobra.Model = None, objective_function_id=None, *args, **kwargs):
        self.reaction_pathway_map = get_reaction_pathway_map(model)
        self.metabolite_pathway_map = get_metabolite_pathway_map(model)
        self.objective_function_id = objective_function_id
        super().__init__(id_or_model=model, *args, **kwargs)

    @classmethod
    def from_sbml(cls, model_path: str, objective_function: str = None):
        model = cobra.io.read_sbml_model(model_path)
        model.objective = objective_function
        return cls(model, objective_function_id=objective_function)

    def __str__(self):
        """
        Print the model summary.
        Returns a string with the model summary.
        -------
        """

        return self.summary()

    @property
    def seeds(self) -> List[Tuple[str, str]]:
        """Identify the seeds of a model.

         Returns:
            A list of tuples, where each tuple contains the ID and compartment of a seed metabolite.
            """

        return self.identify_seeds()

    @property
    def targets(self) -> List[Tuple[str, str]]:
        """Identify the targets of a model.

        Returns:
            A list of tuples, where each tuple contains the ID and compartment of a target metabolite.
            """

        return self.identify_targets()

    # method to identify the seeds of a model
    def identify_seeds(self) -> List[Tuple[str, str]]:
        """Identify the seeds of a model.

        Returns:
            A list of tuples, where each tuple contains the ID and compartment of a seed metabolite.
        """

        total_seeds = []
        total_seeds_ids = []
        for reaction in self.exchanges:
            reactants = reaction.reactants
            for reactant in reactants:
                if reactant.id not in total_seeds_ids and reaction.lower_bound < 0:
                    total_seeds.append((reactant.id, reactant.compartment))
                    total_seeds_ids.append(reactant.id)

        return total_seeds

    def identify_targets(self,
                         objective: str = 'maximize', solver: str = 'cplex') -> List[Tuple[str, str]]:
        """Identify the targets of a model.

        Args:
            objective: The type of objective function. Either 'maximize' or 'minimize'.
            solver: The solver to use. Either 'cplex' or 'glpk'.

        Returns:
            A list of tuples, where each tuple contains the ID and compartment of a target metabolite.
        """

        if not isinstance(self.objective_function_id, str) or not self.objective_function_id:
            objective_function = self.objective
            if not isinstance(objective_function, str) or not objective_function:
                raise ValueError("The objective function must be a string.")

        if not isinstance(objective, str) or not objective:
            raise ValueError("The objective must be a string.")

        if not isinstance(solver, str) or not solver:
            raise ValueError("The solver must be a string.")

        # print(model.summary())
        set_solver(self, solver)

        bio = BioISO(self.objective_function_id, self, objective)
        bio.run(2, False)

        results = bio.get_tree()
        biomass_components = results["M_root_M_root_M_root_product"]["next"]

        targets = []

        for biomass_component in biomass_components:

            biomass_component_role = biomass_components[biomass_component].get("role")
            biomass_component_analysis = biomass_components[biomass_component].get("analysis")
            if biomass_component_role == "Reactant" and not biomass_component_analysis:
                specific_biomass_components = biomass_components[biomass_component].get("next")
                if not specific_biomass_components and biomass_components[biomass_component].get(
                        "identifier") not in targets:
                    targets.append((biomass_components[biomass_component].get("identifier"),
                                    biomass_components[biomass_component].get("compartment")))
                else:
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

    def to_sbml(self, file_name: str, save_path: str, seeds=True, targets=True):
        """
        Write a model to an SBML file.

        Parameters
        ----------
        file_name: str
            The name of the SBML file to write to.
        save_path: str
            The path to save the SBML file to.
        seeds: Boolean
            A list of tuples of seed metabolite IDs and compartments.
        targets: Boolean
            A list of tuples of target metabolite IDs and compartments.
        """
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        if seeds:  # if the user wants to save the seeds.xml
            seeds_file_name = os.path.splitext(file_name)[0] + "_seeds.xml"
            write_metabolites_to_sbml(seeds_file_name, save_path, self.seeds)
        if targets:  # if the user wants to save the targets
            targets_file_name = os.path.splitext(file_name)[0] + "_targets.xml"
            write_metabolites_to_sbml(targets_file_name, save_path, self.targets)

    def create_random_knockout_models(self, num_models: int, min_knockouts: int, max_knockouts: int,
                                      generate_files: bool = False, parent_folder: str = None,
                                      seeds_targets: bool = False, universal_model=None):
        """Create random knockout models.

        Args:
            num_models: The number of models to create.
            min_knockouts: The minimum number of knockouts.
            max_knockouts: The maximum number of knockouts.
            generate_files: Whether to generate SBML files for the models.
            parent_folder: The parent folder to save the models to.
            seeds_targets: Whether to save the seeds and targets of the models.

        Returns:

        """

        # check if the reactions are in the inicial universal model
        universal_model_reactions_ids = set()  # Initialize as an empty set
        kegg_reactions = []

        universal_model_file = None  # Initialize the variable as None

        # If the user wants to generate files and a parent folder was provided
        if universal_model and parent_folder:
            universal_model = cobra.io.read_sbml_model(universal_model)
            universal_model_reactions_ids = set([reaction.id.split('_')[0] for reaction in universal_model.reactions])

            universal_model_file = os.path.join(parent_folder, 'universal_model.xml')
            cobra.io.write_sbml_model(universal_model, universal_model_file)

        # If the user wants to generate files but no parent folder was provided
        elif universal_model:
            universal_model = cobra.io.read_sbml_model(universal_model)
            universal_model_reactions_ids = set([reaction.id.split('_')[0] for reaction in universal_model.reactions])

        # It checks if the reaction ID (without the compartment part) is a subset of the universal_model_reactions_ids
        if universal_model:
            for reaction in self.reactions:
                if {reaction.id.split('_')[0]}.issubset(
                        universal_model_reactions_ids):  # faster when using set instead of list
                    kegg_reactions.append(reaction.id)

        # Generate random knockout models
        knockout_models = []
        knockout_numbers = []

        for i in range(num_models):

            # randomly select the number of knockouts
            num_knockouts = random.randint(min_knockouts, max_knockouts)
            if num_knockouts > len(kegg_reactions):
                num_knockouts = len(kegg_reactions)
                if num_knockouts == 0:
                    raise ValueError(
                        "The number of KEGG reactions is 0. Please check if the KEGG reactions are properly populated.")

            removed_reactions = random.sample(kegg_reactions, num_knockouts)
            new_model = deepcopy(self)  # create a copy of the model
            new_model.remove_reactions(removed_reactions)
            new_model.id = f"{self.id}_knockout_{num_knockouts}"  # set the model ID
            knockout_models.append(new_model)
            knockout_numbers.append(num_knockouts)

        # if the user wants to generate files
        if generate_files and parent_folder:
            for idx, model in enumerate(knockout_models):
                model_folder = os.path.join(parent_folder, f"model_{idx + 1}")
                os.makedirs(model_folder, exist_ok=True)  # exist_ok=True to avoid errors if the folder already exists
                # if the folder was successfully created print a message
                if os.path.exists(model_folder):
                    print(f"Folder {model_folder} created.")

                # write the new model instance to an SBML file
                new_model_instance = Model(model, objective_function_id=self.objective_function_id)
                cobra.io.write_sbml_model(model, os.path.join(model_folder, f"{model.id}.xml"))
                # if the file was successfully created print a message
                if os.path.exists(os.path.join(model_folder, f"{model.id}.xml")):
                    print(f"SBML file {model.id}.xml created.")

                    # create seeds and targets files
                if seeds_targets:
                    new_model_instance.identify_seeds()
                    new_model_instance.identify_targets()
                    new_model_instance.to_sbml(f"{model.id}.xml", save_path=model_folder, seeds=True, targets=True)
                    # if the files were successfully created print a message
                    if os.path.exists(os.path.join(model_folder, f"{model.id}_seeds.xml")):
                        print(f"Seeds file {model.id}_seeds.xml created.")
                    if os.path.exists(os.path.join(model_folder, f"{model.id}_targets.xml")):
                        print(f"Targets file {model.id}_targets.xml created.")

                if universal_model:
                    universal_model_filename = os.path.basename(universal_model_file)
                    dest_universal_model_file = os.path.join(model_folder, universal_model_filename)
                    shutil.copy2(universal_model_file, dest_universal_model_file)
                    if os.path.exists(dest_universal_model_file):
                        print(f"Universal model {universal_model_filename} copied to {model_folder}.")

        return knockout_numbers
