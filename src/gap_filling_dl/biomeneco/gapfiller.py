import os
import time
from functools import partial
from typing import List

import cobra
from bioiso import BioISO, set_solver
from cobra.core import Group
from cobra.io import write_sbml_model
from meneco.query import get_minimal_completion_size
from meneco.meneco import run_meneco
from meneco.meneco import sbml
import copy
from gap_filling_dl.utils import build_temporary_universal_model
from gap_filling_dl.biomeneco.utils import get_compartments_in_model
from gap_filling_dl.biomeneco.utils import identify_dead_end_metabolites
from parallelbar import progress_imap


class GapFiller:
    def __init__(self, model_path, universal_model_path):

        self.universal_model_copy = None
        self.universal_model_compartmentalized = None
        self.dead_ends = None
        self.cloned_model = None
        self.temporary_universal_model = None
        self.original_model = None
        self._universal_model = None
        self.minimal_set = None
        self.minimal_completion = None
        self.gftime = None
        self.cobra_filled_model = None
        self.results_meneco = None
        self.universal_model_path = universal_model_path
        self.model_path = model_path
        self._seeds_path = None
        self._targets_path = None
        self.filled_model = None
        self.objective_function_id = None

    @property
    def seeds_path(self):
        return self._seeds_path

    @seeds_path.setter
    def seeds_path(self, file_path):
        self._seeds_path = file_path

    @property
    def targets_path(self):
        return self._targets_path

    @targets_path.setter
    def targets_path(self, file_path):
        self._targets_path = file_path

    # read the universal model: cobra.io.read_sbml_model(universal_model_path)

    @property
    def universal_model_path(self):
        if self._universal_model_path is None:
            self._universal_model_path = os.path.join(os.path.dirname(self.model_path), 'universal_model.xml')
        return self._universal_model_path

    @universal_model_path.setter
    def universal_model_path(self, value):
        self._universal_model_path = value

    @property
    def universal_model(self):
        if self._universal_model is None:
            self._universal_model = cobra.io.read_sbml_model(self.universal_model_path)
        return self._universal_model

    @universal_model.setter
    def universal_model(self, value):
        self._universal_model = value

    @property
    def original_model(self):
        if self._original_model is None:
            self._original_model = cobra.io.read_sbml_model(self.model_path)
        return self._original_model

    @original_model.setter
    def original_model(self, value):
        self._original_model = value

    def run(self, optimize: bool = False, write_to_file: bool = False, removed_reactions: list = None,
            objective_function_id=None, **kwargs):
        """
        Run the gap-filling algorithm.

        Parameters
        ----------
        optimize
        kwargs:
            Keyword arguments to pass to the gap-filling algorithm.

        """
        if self.original_model is None:
            self.original_model = cobra.io.read_sbml_model(self.model_path)

        if optimize:
            time_start = time.time()
            self.run_meneco(**kwargs)
            time_end = time.time()
            self.gftime = time_end - time_start
            print("Time taken: ", self.gftime)

        else:
            time_start = time.time()

            draftnet = sbml.readSBMLnetwork(self.model_path, 'draft')
            repairnet = sbml.readSBMLnetwork(self.universal_model_path, 'repair')
            seeds = sbml.readSBMLseeds(self.seeds_path)
            targets = sbml.readSBMLtargets(self.targets_path)
            print("Getting minimal completion...")
            # Call the get_minimal_completion_size function this order: draftnet, repairnet, seeds, targets)
            self.minimal_completion = get_minimal_completion_size(draftnet, repairnet, seeds, targets)
            time_end = time.time()
            self.gftime = time_end - time_start

            self.minimal_set = set(
                atom[0] for pred in self.minimal_completion if pred == 'xreaction'
                for atom in self.minimal_completion[pred])

            print("Time taken: ", self.gftime)

        return self.report(write_to_file=write_to_file, removed_reactions=removed_reactions,
                           objective_function_id=objective_function_id, **kwargs)

    def run_meneco(self, enumeration: bool, json_output: bool) -> dict:
        """
        Run the meneco gap-filling algorithm on the model.

        Parameters
        ----------
        enumeration : bool
            If True, the algorithm will enumerate alternative gap-filling solutions.
        json_output : bool
            If True, the algorithm will output the gap-filling solution in JSON format.

        Returns
        -------
        results : dict
            A dictionary containing information about the gap-filling solution.
        """

        self.results_meneco = run_meneco(self.model_path, self.seeds_path, self.targets_path, self.universal_model_path,
                                         enumeration=enumeration, json_output=json_output)

        return self.results_meneco

    @classmethod
    def from_folder(cls, folder_path, temporary_universal_model: bool = False, objective_function_id: str = None,
                    compartments=None, **kwargs):
        """
        Create a GapFiller object from a folder.

        Parameters
        ----------
        objective_function_id
        temporary_universal_model
        folder_path: str
            The path to the folder to create a GapFiller object from.
        """

        model_file = None
        seeds_file = None
        targets_file = None
        universal_model_file = None

        for file in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file)
            print(f"Checking file: {file_path}")
            if file.endswith(".xml"):
                if "model" in file and "universal" not in file and "seeds" not in file and "targets" not in file:
                    model_file = file
                if "seeds" in file:
                    seeds_file = file
                if "targets" in file:
                    targets_file = file
                elif "universal_model" in file and 'temporary' not in file:
                    universal_model_file = file

                print(f"Found file: {file}")

        print(f"model_file: {model_file}")
        print(f"seeds_file: {seeds_file}")
        print(f"targets_file: {targets_file}")
        print(f"universal_model_file: {universal_model_file}")

        # assure that all files are found and raise an error if not
        if not model_file:
            raise FileNotFoundError("No model file found in folder.")
        if not seeds_file:
            raise FileNotFoundError("No seeds file found in folder.")
        if not targets_file:
            raise FileNotFoundError("No targets file found in folder.")
        if not universal_model_file:
            raise FileNotFoundError("No universal model file found in folder.")

        # Read the original model from the model file
        model_path = os.path.join(folder_path, model_file)
        original_model = cobra.io.read_sbml_model(model_path)

        # Create the GapFiller object with the original model
        gap_filler = cls(model_path=model_path, universal_model_path=os.path.join(folder_path, universal_model_file))
        gap_filler.original_model = original_model
        gap_filler.seeds_path = os.path.join(folder_path, seeds_file)
        gap_filler.targets_path = os.path.join(folder_path, targets_file)
        gap_filler.clone_model(compartments=compartments, folder_path=folder_path)
        if temporary_universal_model and gap_filler.temporary_universal_model is None:
            if not objective_function_id:
                raise ValueError(
                    "Objective function ID must be specified for the creation of a temporary universal model.")
            build_temporary_universal_model(gap_filler, folder_path, True)

            if os.path.isfile(os.path.join(folder_path, 'temporary_universal_model.xml')):
                print('Temporary universal model file successfully created.')
                gap_filler.universal_model_path = os.path.join(folder_path, 'temporary_universal_model.xml')
                print('Temporary universal model file path:', gap_filler.universal_model_path)

                gap_filler.temporary_universal_model = cobra.io.read_sbml_model(
                    os.path.join(folder_path, 'temporary_universal_model.xml'))
            else:
                raise FileNotFoundError("Temporary universal model file not found.")
        return gap_filler

    def add_reactions_to_model(self, model, reactions: list[str]):
        """
        Add reactions to a model and run the cobra gap-filling algorithm.

        Parameters
        ----------
        model
        reactions : list
            A list of reactions to add to the model.
        """
        # print('introduced reactions:', reactions)
        # reactions = [reaction.replace('R_', '') if 'R_' in reaction else reaction for reaction in reactions]
        # print('after replacement:', reactions)
        #
        # for r, _ in self.minimal_completion['xreaction']:  # unpack the tuple and ignore the second element
        #     # change the reaction ID to the one used in the universal model
        #     r = r.replace('R_', '')
        #
        #     # search for the reaction in the universal model
        #     reaction = self.universal_model.reactions.get_by_id(r)
        #     # add the reaction to the model
        #     model.add_reactions([reaction])
        #     if reaction in model.reactions:
        #         print('Reaction {} added to model.'.format(r))
        #     else:
        #         print('Reaction {} could not be added to model.'.format(r))
        for reaction in reactions:
            reaction_id = reaction[0]  # Get the reaction ID from the tuple
            # Remove the 'R_' prefix
            clean_reaction_id = reaction_id.replace('R_', '', 1)
            try:
                save_reaction = self.universal_model.reactions.get_by_id(clean_reaction_id)
                model.add_reactions([save_reaction])
            except KeyError:
                print(f"The reaction ID '{clean_reaction_id}' could not be found in the universal model.")

        return model

    def report(self, write_to_file=False, removed_reactions: list = None, objective_function_id: str = None, **kwargs):

        if self.filled_model is None:
            self.filled_model = copy.deepcopy(self.original_model)

        # Include the file names used in the gapfilling without the full path
        report = 'Gap-filling report for {}\n'.format(os.path.basename(self.model_path))
        report += 'Seeds file: {}\n'.format(os.path.basename(self.seeds_path))
        report += 'Targets file: {}\n'.format(os.path.basename(self.targets_path))
        report += 'Universal model file: {}\n'.format(os.path.basename(self.universal_model_path))

        if self.gftime is not None:
            report += "Execution Time: {:.2f} seconds\n".format(self.gftime)

        if self.minimal_completion is not None:
            report += 'This report is based on the minimal completion of the draft network, which is faster. ' \
                      '(This is not the optimized solution.)\n'

            # define the objective function of the model
            if objective_function_id is not None:
                self.filled_model.objective = objective_function_id

            unique_reactions = set(self.minimal_completion['xreaction'])

            # Report the reactions comparison results
            report += '--- Comparison of Reactions ---\n'

            report += 'Added reactions from the gap-filling: {}\n'.format(len(unique_reactions))

            for r in unique_reactions:
                # add the reaction to the model

                report += '- {}\n'.format(r)
                self.add_reactions_to_model(self.filled_model, [r])
                # print the reaction added

                # verify if the reaction is in the model
                if r in self.filled_model.reactions:
                    print('Added reaction: {}'.format(r[0]))

                    # for r in self.minimal_completion['xreaction']:
            #     # add the reaction to the model
            #     self.add_reactions_to_model(self.filled_model, [r])
            #     # print the reaction added
            #     print('Added reaction: {}'.format(r))

            solution = self.filled_model.optimize()

            for r in self.minimal_completion['xreaction']:
                # Remove the 'R_' prefix from the reaction ID
                r_id = r[0].replace('R_', '', 1)

                # get information regarding the reaction added, id, name, lower bound, upper bound, genes and add it
                report += f"\n\nReaction: {r_id}\n"
                report += f"Flux: {solution.fluxes[r_id]}\n"  # Use the solution to access the flux
                report += f"Lower bound: {self.filled_model.reactions.get_by_id(r_id).lower_bound}\n"
                report += f"Upper bound: {self.filled_model.reactions.get_by_id(r_id).upper_bound}\n"
                report += f"Genes: {', '.join([gene.id for gene in self.filled_model.reactions.get_by_id(r_id).genes])}\n"

            # print the fluxes of the reactions
            report += f"Fluxes: {solution.fluxes}\n"

        else:
            report += 'No reactions added from the minimal_completion method of meneco.'

        # if self.results_meneco is not None:
        #     report += "Meneco optimized algorithm results:\n"
        #
        #     if objective_function_id is not None:
        #         self.filled_model.objective = objective_function_id

        #     for model_name, results in self.results_meneco.items():
        #         report += "--- {} ---\n".format(model_name)
        #         report += "Results for {}:\n".format(model_name)
        #         report += "Draft network file: {}\n".format(results["draft_network"])
        #         report += "Seeds file: {}\n".format(results["seeds"])
        #         report += "Targets file: {}\n".format(results["targets"])
        #         report += "Unproducible targets: {}\n".format(results["unproducible_targets"])
        #         report += "Reconstructable targets: {}\n".format(results["reconstructable_targets"])
        #         report += "Essential reactions: {}\n".format(results["essential_reactions"])
        #
        #     if self.results_meneco["One minimal completion"] is not None:
        #         for r in self.results_meneco['One minimal completion']:
        #             # add the reaction to the model
        #             self.add_reactions_to_model(self.filled_model, [r])
        #             # print the reaction added
        #             print('Added reaction: {}'.format(r))
        #
        #         # get information regarding the reaction added, id, name, lower bound, upper bound, genes and add it
        #         report += f"\n\nReaction: {r}\n"
        #         report += f"Flux: {self.filled_model.reactions.get_by_id(r).flux}\n"
        #         report += f"Lower bound: {self.filled_model.reactions.get_by_id(r).lower_bound}\n"
        #         report += f"Upper bound: {self.filled_model.reactions.get_by_id(r).upper_bound}\n"
        #         report += f"Genes: {', '.join([gene.id for gene in self.filled_model.reactions.get_by_id(r).genes])}\n"
        #
        #     solution = self.filled_model.optimize()
        #     # print the fluxes of the reactions
        #     report += f"Fluxes: {solution.fluxes}\n"
        # else:
        #     report += 'No reactions added from the meneco optimized algorithm.'

        # reaçoes removidas antes do gapfilling
        if removed_reactions is not None:
            report += "Previously removed reactions:\n"
            report += "\n".join(removed_reactions)

            # use Cobra package to see fluxes etc of the model
            for reaction in removed_reactions:
                # add reactions to the model
                self.filled_model.add_reactions([self.universal_model.reactions.get_by_id(reaction)])

                # universal_model.reactions.get_by_id('

                report += f"{self.filled_model.reactions.get_by_id(reaction)}\n"
                report += f"Flux: {self.filled_model.reactions.get_by_id(reaction).flux}\n"
                report += f"Lower bound: {self.filled_model.reactions.get_by_id(reaction).lower_bound}\n"
                report += f"Upper bound: {self.filled_model.reactions.get_by_id(reaction).upper_bound}\n"
                report += f"Genes: {', '.join([gene.id for gene in self.filled_model.reactions.get_by_id(reaction).genes])}\n"

        if write_to_file:
            # file_name = report_ + basename.txt
            file_name = os.path.join(os.path.dirname(self.model_path),
                                     "report_" + os.path.basename(self.model_path) + ".txt")
            # check if file exists
            if os.path.isfile(file_name):
                print(f"Warning: {file_name} already exists and will be overwritten.")
            with open(file_name, 'w') as f:
                f.write(report)
            print("Report written to {}".format(file_name))

        return report

    def check_reactions(self, reactions: list):
        """

        Parameters
        ----------
        reactions
        """
        if self.filled_model is None:
            self.filled_model = copy.deepcopy(self.original_model)

        # remove the "R_" prefix from all the reactions if needed
        for reaction in reactions:
            if 'R_' in reaction:
                reactions[reactions.index(reaction)] = reaction.replace('R_', '')

        for reaction in reactions:
            if reaction not in self.universal_model.reactions:
                print(f"Reaction {reaction} not found in the universal model.")
        #     else:
        #         # add reactions to the model
        #         self.filled_model.add_reactions([self.universal_model.reactions.get_by_id(reaction)])
        #         print(f"\n\nReaction: {reaction} has been added to the model.")
        # return self.filled_model

    def clone_model(self, compartments: list = None, **kwargs):
        # Get the compartments in the model
        if not compartments:
            compartments = get_compartments_in_model(self.original_model)

        model_to_use = self.universal_model

        new_reactions = []
        self.universal_model_copy = model_to_use.copy()

        self.universal_model_compartmentalized = cobra.Model()

        self.clone_metabolites(compartments)

        self.clone_reactions(compartments)

        self.clone_groups(compartments)

        self.universal_model_path = os.path.join(kwargs['folder_path'], 'universal_model_compartmentalized.xml')

        write_sbml_model(self.universal_model_compartmentalized, self.universal_model_path)

    def clone_reaction(self, compartments, reactions):
        cloned_reactions = []
        for reaction in reactions:
            for compartment in compartments:
                cloned_reaction = reaction.copy()
                cloned_reaction.id = '__'.join(reaction.id.split("__")[:-1]) + '__' + compartment
                for metabolite in cloned_reaction.metabolites:
                    st = cloned_reaction.metabolites[metabolite]
                    new_metabolite_id_in_compartment = '__'.join(metabolite.id.split("__")[:-1]) + '__' + compartment
                    metabolite_in_compartment = self.universal_model_copy.metabolites.get_by_id(new_metabolite_id_in_compartment)
                    cloned_reaction.add_metabolites({metabolite: -st, metabolite_in_compartment: st})
                cloned_reactions.append(cloned_reaction)
        return cloned_reactions

    def clone_metabolites(self, compartments):
        # print(len(self.universal_model_copy.metabolites))
        # new_metabolites = Parallel(n_jobs=os.cpu_count())(delayed(self.clone_metabolite)(compartments, metabolite) for metabolite in self.universal_model_copy.metabolites)
        list_of_batches = [self.universal_model_copy.metabolites[i:i + 100] for i in range(0, len(self.universal_model_copy.metabolites), 100)]
        new_metabolites = progress_imap(partial(self.clone_metabolite, compartments), list_of_batches)
        new_metabolites = [item for sublist in new_metabolites for item in sublist]
        self.universal_model_copy.add_metabolites(new_metabolites)

    def clone_metabolite(self, compartments, metabolites):
        cloned_metabolites = []
        for metabolite in metabolites:
            for compartment in compartments:
                cloned_metabolite = metabolite.copy()
                cloned_metabolite.id = '__'.join(metabolite.id.split("__")[:-1]) + '__' + compartment
                cloned_metabolite.compartment = compartment
                cloned_metabolites.append(cloned_metabolite)
        return cloned_metabolites

    def get_transport_reactions(self) -> list:
        """
        Get all the transport reactions in the model
        """

        # get all the transport reactions in the model
        transport_reactions = [r for r in self.original_model.reactions if len(r.compartments) > 1]

        # clone the transport reactions
        cloned_transport_reactions = [copy.deepcopy(rxn) for rxn in transport_reactions]

        return cloned_transport_reactions

    def get_dead_end_metabolites(self, model: cobra.Model = None) -> list:
        """
        Identify dead-end metabolites in the model
        """

        # if no model is provided, use the original model and save the result to self.dead_ends
        if model is None:
            model = self.original_model
            self.dead_ends = identify_dead_end_metabolites(model)
            return self.dead_ends
        else:
            # If a model is provided, just return the result without saving it to self.dead_ends
            return identify_dead_end_metabolites(model)

    def identify_dead_end_metabolites_with_bioiso(self, objective_function_id, objective: str = 'maximize',
                                                  solver: str = 'cplex') -> List[str]:
        """
        Identify dead-end metabolites in a model using BioISO.

        Args:
            objective: The type of objective function. Either 'maximize' or 'minimize'.
            solver: The solver to use. Either 'cplex' or 'glpk'.

        Returns:
            A list of dead-end metabolites.

        Parameters
        ----------
        objective_function_id
        """

        self.objective_function_id = objective_function_id

        if not isinstance(self.objective_function_id, str) or not self.objective_function_id:
            objective_function = self.objective_function_id
            if not isinstance(objective_function, str) or not objective_function:
                raise ValueError("The objective function must be a string.")

        if not isinstance(objective, str) or not objective:
            raise ValueError("The objective must be a string.")

        if not isinstance(solver, str) or not solver:
            raise ValueError("The solver must be a string.")

        # Set the solver
        set_solver(self, solver)

        # Run BioISO
        bio = BioISO(self.objective_function_id, self, objective)
        bio.run(2, False)

        # Get the results
        results = bio.get_tree()
        biomass_components = results["M_root_M_root_M_root_product"]["next"]

        dead_ends = []

        for biomass_component in biomass_components:
            biomass_component_role = biomass_components[biomass_component].get("role")
            biomass_component_analysis = biomass_components[biomass_component].get("analysis")
            if biomass_component_role == "Reactant" and not biomass_component_analysis:
                dead_ends.append(biomass_components[biomass_component].get("identifier"))

        return dead_ends

    def fill_dead_ends(self):
        """
        Fill the dead-end metabolites in the model using transport reactions.
        """
        if self.dead_ends is None:
            raise ValueError(
                "Dead-end metabolites have not been identified. Please run the dead-end identification step first.")
        if self.cloned_model is None:
            raise ValueError("The model has not been cloned. Please run the cloning step first.")

        # Filter dead-end metabolites that are already present in the model
        filtered_dead_ends = [metabolite for metabolite in self.dead_ends if
                              metabolite not in self.cloned_model.metabolites]

        transport_reactions = []
        for metabolite_id in filtered_dead_ends:
            transport_reaction = self.create_and_add_transport_reaction(metabolite_id)
            transport_reactions.append(transport_reaction.id)

        print(f"Added transport reactions to move {len(transport_reactions)} dead-end metabolites.")

        return transport_reactions

    def get_plausible_transports_for_dead_ends(self):
        """
        Create a dictionary of plausible transport reactions for each dead-end metabolite.
        """
        if self.dead_ends is None:
            raise ValueError(
                "Dead-end metabolites have not been identified. Please run the dead-end identification step first.")

        # A dictionary to store the plausible transport reactions for each dead-end metabolite
        plausible_transports = {}

        # Identify transport reactions in the cloned model
        for reaction in self.cloned_model.reactions:
            # Transport reactions involve metabolites in different compartments
            if len(reaction.compartments) > 1:
                for metabolite in reaction.metabolites:
                    # If the metabolite is a dead-end metabolite, add the reaction to its plausible transports
                    if metabolite in self.dead_ends:
                        if metabolite not in plausible_transports:
                            plausible_transports[metabolite] = []
                        plausible_transports[metabolite].append(reaction)

        return plausible_transports

    def create_and_add_transport_reaction(self, metabolite_id):
        """
        Create a new transport reaction for the given metabolite and add it to the cloned model.
        """
        # Create a new transport reaction
        reaction_id = "TR_" + metabolite_id
        transport_reaction = cobra.Reaction(reaction_id)
        transport_reaction.name = "Transport reaction for metabolite {}".format(metabolite_id)

        # Add the metabolite to transport from its current compartment to a different compartment
        metabolite = self.cloned_model.metabolites.get_by_id(metabolite_id)
        transport_reaction.add_metabolites({metabolite: -1.0})
        transport_reaction.add_metabolites({metabolite_id + "_transport": 1.0})

        # If the transport reaction is plausible, add it to the model
        if self.is_plausible_transport(transport_reaction):
            transport_reaction.lower_bound = -1000.0  # Allow negative flux (transport out of the compartment)
            transport_reaction.upper_bound = 1000.0  # Allow positive flux (transport into the compartment)

            # Add the reaction to the cloned model
            self.cloned_model.add_reaction(transport_reaction)

        return transport_reaction

    def is_plausible_transport(self, transport_reaction):
        """
        Check if the given transport reaction is plausible.
        """
        # Check if the reaction is balanced
        for metabolite in transport_reaction.metabolites:
            if abs(transport_reaction.get_coefficient(metabolite)) != 1.0:
                return False

        # Check if the reaction is reversible, there are trnsport reactions that are irreversible tho
        if transport_reaction.lower_bound < 0.0 or transport_reaction.upper_bound < 0.0:
            return False

        # Check if the reaction is thermodynamically feasible

    def clone_reactions(self, compartments):
        """
        Clone reactions in the model.
        Parameters
        ----------
        compartments

        Returns
        -------

        """
        # results = Parallel(n_jobs=4)(delayed(self.clone_reaction)(reaction, compartment) for reaction in self.universal_model_copy.reactions for compartment in compartments)
        list_of_batches = [self.universal_model_copy.reactions[i:i + 100] for i in
                           range(0, len(self.universal_model_copy.reactions), 100)]
        reactions = progress_imap(partial(self.clone_reaction, compartments), list_of_batches)
        reactions = [item for sublist in reactions for item in sublist]
        self.universal_model_compartmentalized.add_reactions(reactions)

    def clone_groups(self, compartments):
        """
        Clone groups in the model.
        Returns
        -------

        """
        all_reactions_ids = [reaction.id for reaction in self.universal_model_compartmentalized.reactions]
        new_groups = []
        for group in self.universal_model_copy.groups:
            new_members = []
            group_copy = Group(id=group.id, name=group.name)
            for compartment in compartments:
                for member in group.members:
                    reaction_id = '__'.join(member.id.split("__")[:-1]) + '__' + compartment
                    if reaction_id in all_reactions_ids:
                        new_members.append(self.universal_model_compartmentalized.reactions.get_by_id(reaction_id))
                    else:
                        print(f"Reaction {reaction_id} not found in the model")
            group_copy.add_members(new_members)
            new_groups.append(group_copy)
        self.universal_model_compartmentalized.add_groups(new_groups)
