import copy
import os
import random
import re
import shutil
from copy import deepcopy
from functools import partial
from typing import Tuple, List

import cobra
from bioiso import BioISO
from bioiso import set_solver
from cobra import Reaction
from pandas import DataFrame
from parallelbar import progress_imap

from .utils import get_metabolite_pathway_map, get_reaction_pathway_map

try:
    from write_sbml import write_metabolites_to_sbml
except:
    from gap_filling_dl.write_sbml.metabolites_to_sbml import write_metabolites_to_sbml


class Model(cobra.Model):
    def __init__(self, model: cobra.Model = None, objective_function_id=None, *args, **kwargs):
        """
        Initialize a model.
        Parameters
        ----------
        model : cobra.Model
            A cobra model.
                objective_function_id : str, optional
            The ID of the objective function to be used in this model.
        *args
            Variable length argument list.
        **kwargs
            Arbitrary keyword arguments.

        Attributes
        ----------
        e_res_precursors : NoneType
            Placeholder for essential resource precursors in the model.
        precursors_reactions : NoneType
            Placeholder for precursor reactions in the model.
        reaction_pathway_map : dict
            Mapping of reactions to their respective pathways in the model.
        metabolite_pathway_map : dict
            Mapping of metabolites to their respective pathways in the model.
        pathway_metabolites_map : dict
            Mapping of pathways to their respective metabolites in the model.
        objective_function_id : str
            The ID of the objective function in the model.
        pre_precursors : NoneType
            Placeholder for pre-precursors in the model.
        bio_precursors : NoneType
            Placeholder for biological precursors in the model.
        _seeds : NoneType
            Placeholder for seeds in the model.
        _targets : NoneType
            Placeholder for targets in the model.
        model_old : list
            Placeholder for storing the original state of the model.
        """
        self.e_res_precursors = None
        self.precursors_reactions = None
        if model and model.reactions and model.groups:
            self.reaction_pathway_map = get_reaction_pathway_map(model)
        if model and model.metabolites and model.groups:
            self.metabolite_pathway_map = get_metabolite_pathway_map(model)
            self.pathway_metabolites_map = {}
            for metabolite, pathways in self.metabolite_pathway_map.items():
                for pathway in pathways:
                    if pathway not in self.pathway_metabolites_map:
                        self.pathway_metabolites_map[pathway] = []
                    self.pathway_metabolites_map[pathway].append(metabolite)
        self.objective_function_id = objective_function_id
        self.pre_precursors = None
        self.bio_precursors = None
        self._seeds = None
        self._targets = None
        self.model_old = []
        super().__init__(id_or_model=model, *args, **kwargs)

    @classmethod
    def from_sbml(cls, model_path: str, objective_function: str = None):
        """
        Class method to create a Model instance from an SBML file.

        Parameters
        ----------
        model_path : str
            Path to the SBML file.
        objective_function : str, optional
            The ID of the objective function to be used in this model. If not provided, the model's current objective function will be used.

        Returns
        -------
        Model
            A Model instance initialized with the data from the SBML file and the specified objective function.
        """
        model = cobra.io.read_sbml_model(model_path)
        model.objective = objective_function
        return cls(model, objective_function_id=objective_function)

    def __str__(self):
        """
        Print the model summary.

        Returns
        -------
        str
            The model summary.
        """

        return self.summary()

    @property
    def seeds(self) -> List[Tuple[str, str]]:
        """
        Get the seeds of the model.

        This property returns the seeds of the model, which are stored in the `_seeds` attribute. The seeds are represented
        as a list of tuples, where each tuple contains the ID and compartment of a seed metabolite.

        Returns
        -------
        List[Tuple[str, str]]
            The seeds of the model. Each seed is represented as a tuple containing the ID and compartment of a seed metabolite.
        """
        return self._seeds

    @seeds.setter
    def seeds(self, seeds: List[Tuple[str, str]]):
        """
        Set the seeds of the model.

        This method sets the seeds of the model, which are stored in the `_seeds` attribute. The seeds are represented
        as a list of tuples, where each tuple contains the ID and compartment of a seed metabolite.

        Parameters
        ----------
        seeds : List[Tuple[str, str]]
            The seeds to set for the model. Each seed is represented as a tuple containing the ID and compartment of a seed metabolite.
        """
        self._seeds = seeds

    @property
    def targets(self) -> List[Tuple[str, str]]:
        """
        Get the targets of the model.

        This property returns the targets of the model, which are stored in the `_targets` attribute. The targets are represented
        as a list of tuples, where each tuple contains the ID and compartment of a target metabolite.

        Returns
        -------
        List[Tuple[str, str]]
            The targets of the model. Each target is represented as a tuple containing the ID and compartment of a target metabolite.
        """
        return self._targets

    @targets.setter
    def targets(self, targets: List[Tuple[str, str]]):
        """
        Set the targets of the model.

        This method sets the targets of the model, which are stored in the `_targets` attribute. The targets are represented
        as a list of tuples, where each tuple contains the ID and compartment of a target metabolite.

        Parameters
        ----------
        targets : List[Tuple[str, str]]
            The targets to set for the model. Each target is represented as a tuple containing the ID and compartment of a target metabolite.
        """
        self._targets = targets

    # method to identify the seeds of a model
    def identify_seeds(self) -> List[Tuple[str, str]]:
        """
        Identify the seeds of the model.

        This method identifies the seeds of the model by iterating over the boundary reactions and checking the reactants.
        A seed is defined as a reactant in a boundary reaction with a lower bound less than 0. The seeds are stored in the `_seeds` attribute.

        Returns
        -------
        List[Tuple[str, str]]
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
        self.seeds = total_seeds
        return total_seeds

    def add_seeds(self, new_seeds: List[Tuple[str, str]]):
        """
        Add new seeds to the model.

        This method adds new seeds to the model by extending the `_seeds` attribute with the provided list of seeds.

        Parameters
        ----------
        new_seeds : List[Tuple[str, str]]
            A list of tuples, where each tuple contains the ID and compartment of a seed metabolite to be added to the model.
        """
        self.seeds.extend(new_seeds)

    def identify_targets(self, objective: str = 'maximize', solver: str = 'cplex') -> List[Tuple[str, str]]:
        """
        Identify the targets of a metabolic model.

        This method uses the BioISO algorithm to identify the targets of the model. A target is a metabolite that is produced by the model's reactions. The method takes two optional arguments: `objective` and `solver`. The `objective` argument specifies the type of objective function to use in the BioISO algorithm, either 'maximize' or 'minimize'. The `solver` argument specifies the solver to use in the BioISO algorithm, either 'cplex' or 'glpk'.

        Parameters
        ----------
        objective : str, optional
            The type of objective function to use in the BioISO algorithm. Either 'maximize' or 'minimize'. Default is 'maximize'.
        solver : str, optional
            The solver to use in the BioISO algorithm. Either 'cplex' or 'glpk'. Default is 'cplex'.

        Returns
        -------
        List[Tuple[str, str]]
        A list of tuples, where each tuple contains the ID and compartment of a target metabolite.
        """
        self._validate_inputs(objective, solver)
        set_solver(self, solver)

        bio = BioISO(self.objective_function_id, self, objective, time_out=None)
        bio.run(2, False)

        results = bio.get_tree()
        biomass_components = results["M_root_M_root_M_root_product"]["next"]

        targets = self._extract_targets(biomass_components)

        not_produced = self.test_e_precursors(to_ignore=[e[0] for e in targets])

        for not_produced in not_produced:
            targets.add((not_produced, self.metabolites.get_by_id(not_produced).compartment))
        self.targets = list(targets)
        return list(self.targets)

    def _validate_inputs(self, objective, solver):
        """
            Validate the inputs to the identify_targets method.

    This method checks the validity of the inputs to the identify_targets method. It checks that the objective_function_id attribute,
    the objective argument, and the solver argument are all strings and are not empty. If any of these conditions are not met,
    it raises a ValueError.

    Parameters
    ----------
    objective : str
        The type of objective function to use in the BioISO algorithm. Either 'maximize' or 'minimize'.
    solver : str
        The solver to use in the BioISO algorithm. Either 'cplex' or 'glpk'.

    Raises
    ------
    ValueError
        If the objective_function_id attribute, the objective argument, or the solver argument is not a string or is empty.
        """
        if not isinstance(self.objective_function_id, str) or not self.objective_function_id:
            objective_function = self.objective
            if not isinstance(objective_function, str) or not objective_function:
                raise ValueError("The objective function must be a string.")

        if not isinstance(objective, str) or not objective:
            raise ValueError("The objective must be a string.")

        if not isinstance(solver, str) or not solver:
            raise ValueError("The solver must be a string.")

    def _extract_targets(self, biomass_components):
        """
            Extract the targets from the biomass components.

            This method iterates over the biomass components and identifies the targets. A target is defined as a reactant
            in a biomass component that does not have an associated analysis. The targets are stored in a set to avoid duplicates.

            Parameters
            ----------
            biomass_components : dict
                A dictionary containing the biomass components. Each key is a biomass component and the value is a dictionary
                containing information about the component, such as its role and analysis.

            Returns
            -------
            set
                A set of tuples, where each tuple contains the ID and compartment of a target metabolite.
        """
        targets = set()

        for biomass_component in biomass_components:
            biomass_component_role = biomass_components[biomass_component].get("role")
            biomass_component_analysis = biomass_components[biomass_component].get("analysis")
            if biomass_component_role == "Reactant" and not biomass_component_analysis:
                specific_biomass_components = biomass_components[biomass_component].get("next")
                if not specific_biomass_components and biomass_components[biomass_component].get(
                        "identifier") not in targets:
                    targets.add((biomass_components[biomass_component].get("identifier"),
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
                                targets.add((target_id, compartment))
        return targets

    def test_e_precursors(self, to_ignore=[]):
        """
        Test the essential precursors of the biomass.

    This method tests any precursor of the biomass by evaluating the flux of each precursor. If the flux is zero,
    the precursor is considered essential. The method takes an optional argument `to_ignore` which is a list of
    precursor IDs to ignore during the test.

    Parameters
    ----------
    to_ignore : list, optional
        A list of precursor IDs to ignore during the test. Default is an empty list.

    Returns
    -------
    Index
        An index object (similar to a list) of the IDs of the essential precursors of the biomass.

        """
        if self.bio_precursors is None:
            self.get_bio_precursors()
        if self.pre_precursors is None:
            self.get_pre_precursors()

        res = progress_imap(partial(self.evaluate_preprecursors, to_ignore), self.bio_precursors)
        e_precursors_res = {"Flux": []}
        meta_val = []
        for e_precursor_res, meta in res:
            e_precursors_res["Flux"].extend(e_precursor_res["Flux"])
            meta_val.extend(meta)

        self.e_res_precursors = DataFrame(e_precursors_res, index=meta_val)

        return self.e_res_precursors.loc[self.e_res_precursors["Flux"] == 0].index

    def evaluate_preprecursors(self, to_ignore, precursor):
        """
        Evaluate the precursors of the precursors (pre-precursors) of a given metabolite.

        This method checks if a given precursor is a pre-precursor and not in the list of metabolites to ignore. If it is,
        it creates a temporary copy of the model, evaluates the flux of the pre-precursor in the temporary model, and
        appends the flux value to the results. The method returns the results and the IDs of the evaluated pre-precursors.

        Parameters
        ----------
        to_ignore : list
            A list of metabolite IDs to ignore during the evaluation.
        precursor : cobra.Metabolite
            The metabolite for which to evaluate the pre-precursors.

        Returns
        -------
        dict
            A dictionary with the key "Flux" and a list of flux values as the value.
        list
            A list of the IDs of the evaluated pre-precursors.

        """
        e_precursors_res = {"Flux": []}
        meta_val = []
        if precursor.id not in to_ignore and precursor.id in self.pre_precursors:
            # if precursor.id not in self.pre_precursors:
            #     self.pre_precursors[precursor.id] = []
            for pre_precursor in self.pre_precursors[precursor.id]:
                if pre_precursor.id not in to_ignore:
                    temp_model = copy.deepcopy(self)
                    e_precursors_res["Flux"].append(self.evaluate_precursor(temp_model, pre_precursor))
                    meta_val.append(pre_precursor.id)
        return e_precursors_res, meta_val

    def evaluate_precursor(self, temp_model, current_precursor):
        """
        Evaluate the flux of a given precursor in a given model.

        This method creates a demand reaction for the given precursor in the given model, sets the model's objective to
        the demand reaction, and then optimizes the model to calculate the flux. The method returns the flux value.

        Parameters
        ----------
        temp_model : Model
            The model in which to evaluate the precursor.
        current_precursor : cobra.Metabolite
            The precursor to evaluate.

        Returns
        -------
        float
        The flux value of the precursor in the model.

        """
        reaction = temp_model.create_demand(current_precursor.id)
        temp_model.objective = reaction.id
        val = temp_model.slim_optimize()
        return val

    def to_sbml(self, file_name: str, save_path: str, seeds=True, targets=True):
        """
        Save the model, seeds, and targets to SBML files.

        This method saves the current model, its seeds, and targets to separate SBML files. The files are saved in the specified directory. If the directory does not exist, it is created. The names of the seeds and targets files are derived from the name of the model file by appending '_seeds' and '_targets' respectively.

        Parameters
        ----------
        file_name : str
            The name of the SBML file to which the model will be saved.
        save_path : str
            The path to the directory where the SBML files will be saved.
        seeds : bool, optional
            Whether to save the seeds of the model to an SBML file. Default is True.
        targets : bool, optional
            Whether to save the targets of the model to an SBML file. Default is True.

        Returns
        -------
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
        """
        Create random knockout models.

        This method creates a specified number of random knockout models by randomly removing a specified number of reactions from the model. The method can also generate SBML files for the knockout models and save them to a specified directory. The method can also include the seeds and targets of the model in the SBML files.

        Parameters
        ----------
        num_models : int
            The number of knockout models to create.
        min_knockouts : int
            The minimum number of reactions to remove from the model to create a knockout model.
        max_knockouts : int
            The maximum number of reactions to remove from the model to create a knockout model.
        generate_files : bool, optional
            Whether to generate SBML files for the knockout models. Default is False.
        parent_folder : str, optional
            The directory where the SBML files will be saved. Required if generate_files is True. Default is None.
        seeds_targets : bool, optional
            Whether to include the seeds and targets of the model in the SBML files. Default is False.
        universal_model : str, optional
            The path to the universal model. Default is None.

        Returns
        -------
        list
            A list of the number of knockouts in each generated knockout model.
        """
        # Initialize the universal model reactions IDs and the KEGG reactions
        universal_model_reactions_ids, kegg_reactions, universal_model_file = self._initialize_universal_model_reactions_ids(
            universal_model, parent_folder)

        # Generate random knockout models
        knockout_models, knockout_numbers, removed_reactions_dict = self._generate_random_knockout_models(
            num_models, min_knockouts, max_knockouts, kegg_reactions)

        # Generate files if requested
        if generate_files and parent_folder:
            self._generate_files(knockout_models, parent_folder, seeds_targets, universal_model_file)

        # Print the reactions removed from each model
        self._print_removed_reactions(removed_reactions_dict)

        return knockout_numbers

    def _initialize_universal_model_reactions_ids(self, universal_model, parent_folder):
        """
        Initialize the universal model reactions IDs and the KEGG reactions.

        This method checks if the reactions are in the initial universal model. If the user wants to generate files and a parent folder was provided, it writes the universal model to an SBML file.

        Parameters
        ----------
        universal_model : str
            The path to the universal model.
        parent_folder : str
            The directory where the SBML files will be saved.

        Returns
        -------
        set
            A set of the universal model reactions IDs.
        list
            A list of the KEGG reactions.
        str
            The path to the universal model file.
        """
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

        return universal_model_reactions_ids, kegg_reactions, universal_model_file

    def _generate_random_knockout_models(self, num_models, min_knockouts, max_knockouts, kegg_reactions):
        """
        Generate random knockout models.

        This method generates a specified number of random knockout models by randomly removing a specified number of reactions from the model.

        Parameters
        ----------
        num_models : int
            The number of knockout models to create.
        min_knockouts : int
            The minimum number of reactions to remove from the model to create a knockout model.
        max_knockouts : int
            The maximum number of reactions to remove from the model to create a knockout model.
        kegg_reactions : list
            A list of the KEGG reactions.

        Returns
        -------
        list
            A list of the generated knockout models.
        list
            A list of the number of knockouts in each generated knockout model.
        dict
            A dictionary of the reactions removed from each model.
        """
        knockout_models = []
        knockout_numbers = []
        removed_reactions_dict = {}
        max_attempts = num_models * 10
        generated_models = 0
        attempts = 0

        while generated_models < num_models:
            attempts += 1
            if attempts > max_attempts:
                raise Exception('Unable to generate the desired number of models')
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
                continue # Skip the current iteration and do not add the model to the list

            knockout_models.append(new_model)
            knockout_numbers.append(num_knockouts)
            removed_reactions_dict[generated_models] = removed_reactions
            generated_models += 1

        return knockout_models, knockout_numbers, removed_reactions_dict

    def _generate_files(self, knockout_models, parent_folder, seeds_targets, universal_model_file):
        for idx, model in enumerate(knockout_models):
            model_folder = os.path.join(parent_folder, f"model_{idx + 1}")
            os.makedirs(model_folder, exist_ok=True)

            # write the new model instance to an SBML file
            new_model_instance = Model(model, objective_function_id=self.objective_function_id)
            cobra.io.write_sbml_model(model, os.path.join(model_folder, f"{model.id}.xml"))

            # create seeds and targets files
            if seeds_targets:
                new_model_instance.identify_seeds()
                new_model_instance.identify_targets()
                new_model_instance.to_sbml(f"{model.id}.xml", save_path=model_folder, seeds=True, targets=True)

            if universal_model_file:
                universal_model_filename = os.path.basename(universal_model_file)
                dest_universal_model_file = os.path.join(model_folder, universal_model_filename)
                shutil.copy2(universal_model_file, dest_universal_model_file)

    def _print_removed_reactions(self, removed_reactions_dict):
        """
        Print the reactions removed from each model.

        This method prints the reactions removed from each model.

        Parameters
        ----------
        removed_reactions_dict : dict
            A dictionary of the reactions removed from each model.
        """
        print("Reactions removed from each model:")
        for key, value in removed_reactions_dict.items():
            print(f"Model {key + 1}: {value}")

    def get_bio_precursors(self):
        """
        Retrieve the precursors of the biomass reaction.

        This method identifies the reactants of the objective function of the model, which are considered as the precursors of the biomass reaction. The identified precursors are stored in the `bio_precursors` attribute of the model.

        Returns
        -------
        list
            A list of cobra.Metabolite objects representing the precursors of the biomass reaction.
        """
        self.bio_precursors = self.get_reactants(self.objective_function_id)
        return self.bio_precursors

    def get_pre_precursors(self):
        """
        Retrieve the precursors of the biomass precursors.

        This method identifies the reactants of the reactions in which the biomass precursors participate. These reactants are considered as the precursors of the biomass precursors. The identified precursors are stored in the `pre_precursors` attribute of the model.

        Returns
        -------
        dict
            A dictionary where the keys are the IDs of the biomass precursors and the values are lists of cobra.Metabolite objects representing the precursors of the corresponding biomass precursor.
        """
        try:
            self.pre_precursors = {}
            self.precursors_reactions = {}
            if self.bio_precursors is None:
                self.get_bio_precursors()
            for precursor in self.bio_precursors:
                reactions = precursor.reactions
                if len(reactions) > 2 and re.match(r'^[C]\d{5}__.*', precursor.id):
                    self.pre_precursors[precursor.id] = []
                    print(f"{precursor.id} is a product of more than one reaction")
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
        """Retrieve a reaction from the model by its ID.

        This method attempts to retrieve a reaction from the model by its ID. If the reaction is not found, it returns False.

        Parameters
        ----------
        reaction : str
            The ID of the reaction to retrieve.

        Returns
        -------
        cobra.Reaction or bool
        The reaction object if found, otherwise False.
        """
        try:
            return self.reactions.get_by_id(reaction)
        except Exception as e:
            return False

    def get_products(self, reaction):
        """
        Retrieve the products of a reaction.

        This method retrieves the products of a given reaction in the model.

        Parameters
        ----------
        reaction : str
            The ID of the reaction.

        Returns
        -------
        list
            A list of cobra.Metabolite objects representing the products of the reaction.
        """
        return self.get_reaction(reaction).products

    def create_trnas_reactions(self, protein_id="e-Protein"):
        """
        Create tRNA reactions for protein synthesis.

        This method creates tRNA reactions for protein synthesis in the model. It identifies the products of the protein synthesis reaction that are not water or protein, and creates a sink reaction for each of them.

        Parameters
        ----------
        protein_id : str, optional
            The ID of the protein synthesis reaction. Default is "e-Protein".

        Returns
        -------
        None
    """
        print()
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
        """
        Create a demand reaction for a metabolite.

        This method creates a demand reaction for a given metabolite in the model. If a demand reaction for the metabolite already exists, it does not create a new one.

        Parameters
        ----------
        metabolite_id : str
            The ID of the metabolite for which to create a demand reaction.

        Returns
        -------
        cobra.Reaction
            The demand reaction for the metabolite.
        """

        if not self.get_reaction("DM_" + metabolite_id):
            reaction_name = "DM_" + metabolite_id
            self.create_reaction(reaction_name).add_metabolites({self.get_metabolite(metabolite_id): -1})
            self.get_reaction(reaction_name).bounds = (0, 10000)
            return self.get_reaction(reaction_name)

    def create_sink(self, metabolite_id, bounds=(-10000, 10000)):
        """
        Create a sink reaction for a metabolite.

        This method creates a sink reaction for a given metabolite in the model. If a sink reaction for the metabolite already exists, it does not create a new one.

        Parameters
        ----------
        metabolite_id : str
            The ID of the metabolite for which to create a sink reaction.
        bounds : tuple, optional
            The lower and upper bounds for the sink reaction. Default is (-10000, 10000).

        Returns
        -------
        cobra.Reaction
            The sink reaction for the metabolite.
        """
        if not self.get_reaction("Sk_" + metabolite_id):
            reaction_name = "Sk_" + metabolite_id
            self.create_reaction(reaction_name).add_metabolites({self.get_metabolite(metabolite_id): -1})
            self.get_reaction(reaction_name).bounds = bounds
            return self.get_reaction(reaction_name)

    def create_reaction(self, reaction):
        """
        Create a reaction in the model.

        This method creates a reaction in the model with the given ID.

        Parameters
        ----------
        reaction : str
            The ID of the reaction to create.

        Returns
        -------
        cobra.Reaction
            The created reaction.
        """
        self.add_reactions([cobra.Reaction(reaction)])
        return self.get_reaction(reaction)

    def get_metabolite(self, metabolite):
        """
        Retrieve a metabolite from the model by its ID.

        This method retrieves a metabolite from the model by its ID.

        Parameters
        ----------
        metabolite : str
            The ID of the metabolite to retrieve.

        Returns
        -------
        cobra.Metabolite
            The metabolite object.
        """
        return self.metabolites.get_by_id(metabolite)

    def get_reactants(self, reaction):
        """
        Retrieve the reactants of a reaction.

        This method retrieves the reactants of a given reaction in the model.

        Parameters
        ----------
        reaction : str
            The ID of the reaction.

        Returns
        -------
        list
            A list of cobra.Metabolite objects representing the reactants of the reaction.
        """
        return self.get_reaction(reaction).reactants
