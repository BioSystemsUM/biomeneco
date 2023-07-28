import os
import random
import shutil
from copy import deepcopy
from typing import Tuple, List

import cobra
from bioiso import BioISO
from bioiso import set_solver
from cobra import Reaction

from .utils import get_metabolite_pathway_map, get_reaction_pathway_map
from write_sbml import write_metabolites_to_sbml


class Model(cobra.Model):
    def __init__(self, model: cobra.Model = None, objective_function_id=None, *args, **kwargs):
        self.precursors_reactions = None
        self.reaction_pathway_map = get_reaction_pathway_map(model)
        self.metabolite_pathway_map = get_metabolite_pathway_map(model)
        self.objective_function_id = objective_function_id
        self.pre_precursors = None
        self.bio_precursors = None
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
        for reaction in self.boundary:
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
        Save the model to SBML.

        Args:
            file_name: The name of the SBML file.
            save_path: The path to the folder where the SBML file will be saved.
            seeds: Whether to save the seeds of the model.
            targets: Whether to save the targets of the model.

        Returns:
            None
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
            num_models: The number of knockout models to create.
            min_knockouts: The minimum number of knockouts to create.
            max_knockouts: The maximum number of knockouts to create.
            generate_files: Whether to generate SBML files for the knockout models.
            parent_folder: The parent folder to save the SBML files to.
            seeds_targets: Whether to include the seeds and targets in the SBML files.
            universal_model: The path to the universal model.

        Returns:
            A list of random knockout models.

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
        removed_reactions_dict = {}

        generated_models = 0

        while generated_models < num_models:
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

            # Check if seeds and targets are not empty
            new_model_instance = Model(new_model, objective_function_id=self.objective_function_id)
            new_model_instance.identify_seeds()
            new_model_instance.identify_targets()
            if not new_model_instance.seeds or not new_model_instance.targets:
                continue  # Skip the current iteration and do not add the model to the list

            knockout_models.append(new_model)
            knockout_numbers.append(num_knockouts)
            removed_reactions_dict[generated_models] = removed_reactions
            generated_models += 1

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

        print("Reactions removed from each model:")
        for key, value in removed_reactions_dict.items():
            print(f"Model {key + 1}: {value}")

        return knockout_numbers

    def get_bio_precursors(self):
        """
        This function returns the precursors of the biomass reaction
        Returns
        -------

        """
        self.bio_precursors = self.get_reactants(self.objective_function_id)
        return self.bio_precursors

    def get_pre_precursors(self):

        """ This function returns the precursors of the biomass precursors
        retrieved earlier in the get_bio_precursors() function"""
        try:
            self.pre_precursors = {}
            self.precursors_reactions = {}
            if self.bio_precursors is None:
                self.get_bio_precursors()
            for precursor in self.bio_precursors:
                reactions = precursor.reactions
                if len(reactions) > 2:
                    print("ATP or other metabolite. Without unique precursor")
                    self.pre_precursors[precursor.id] = []
                else:
                    reaction = Reaction()
                    for r in reactions:
                        if r.id != self.objective_function_id:
                            reaction = r
                    if reaction.id:
                        self.pre_precursors[precursor.id] = self.get_reactants(reaction.id)
                        self.precursors_reactions[precursor.name] = (reaction.id, precursor.id)
        except Exception as e:
            print(e)
        return self.pre_precursors

    def get_reaction(self, reaction):
        try:
            return self.reactions.get_by_id(reaction)
        except Exception as e:
            print(e)
            print(reaction + " not found")

    def get_products(self, reaction):
        return self.get_reaction(reaction).products

    def create_trnas_reactions(self, protein_id="e-Protein"):
        """
        This function creates the tRNA reactions for the protein synthesis
        """
        try:
            if self.pre_precursors is None:
                self.get_pre_precursors()

            trnas = self.get_products(self.precursors_reactions[protein_id][0])

            for trna in trnas:
                if "H2O" not in trna.name and "e-Protein" not in trna.name:
                    self.create_sink(trna.id)
        except Exception as e:
            print("No protein found in the model")
            print(e)

    def create_demand(self, metabolite_id):

        """ This function creates a demand reaction"""

        reaction_name = "DM_" + metabolite_id
        self.create_reaction(reaction_name).add_metabolites({self.get_metabolite(metabolite_id): -1})
        self.get_reaction(reaction_name).bounds = (0, 10000)
        return self.get_reaction(reaction_name)

    def create_sink(self, metabolite_id, bounds=(-10000, 10000)):

        """ This function creates a sink reaction """

        reaction_name = "Sk_" + metabolite_id
        self.create_reaction(reaction_name).add_metabolites({self.get_metabolite(metabolite_id): -1})
        self.get_reaction(reaction_name).bounds = bounds
        return self.get_reaction(reaction_name)

    def create_reaction(self, reaction):
        self.add_reactions([cobra.Reaction(reaction)])
        return self.get_reaction(reaction)

    def get_metabolite(self, metabolite):
        return self.metabolites.get_by_id(metabolite)

    def get_reactants(self, reaction):
        return self.get_reaction(reaction).reactants
