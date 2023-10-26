import copy
import itertools
import json
import os
import platform
import shutil
import time
from functools import partial
from os.path import join, dirname
from typing import List

import cobra
import parallelbar
from bioiso import BioISO, set_solver
from clyngor.as_pyasp import TermSet, Atom
from cobra import Reaction, Metabolite
from cobra.core import Group
from cobra.io import write_sbml_model, read_sbml_model
from meneco.meneco import run_meneco, extract_xreactions
from meneco.meneco import sbml
from meneco.query import get_minimal_completion_size
from parallelbar import progress_imap
from parallelbar.parallelbar import ProgressBar, _update_error_bar
from typing_extensions import deprecated
from write_sbml import write_metabolites_to_sbml

from gap_filling_dl.kegg_api import create_related_pathways_map, get_related_pathways
from .model import Model
from .utils import get_compartments_in_model
from .utils import identify_dead_end_metabolites


### The following code is to set up the progress bar. The parallelbar package does not allow changing the description is the bar, so I had to change the code.
parallelbar.parallelbar.PROGRESS = "DONE"
def _my_process_status(bar_size, worker_queue):
    bar = ProgressBar(total=bar_size, desc=parallelbar.parallelbar.PROGRESS)
    error_bar_parameters = dict(total=bar_size, position=1, desc='ERROR', colour='red')
    error_bar = {}
    error_bar_n = 0
    while True:
        flag, upd_value = worker_queue.get()
        if flag is None:
            if error_bar:
                error_bar_n = error_bar['bar'].n
                error_bar['bar'].close()
            if bar.n < bar_size and upd_value != -1:
                bar.update(bar_size - bar.n - error_bar_n)
            bar.close()
            break
        if flag:
            _update_error_bar(error_bar, error_bar_parameters)
        else:
            bar.update(upd_value)


parallelbar.parallelbar._process_status = _my_process_status



if platform.system() == 'Linux':
    UTILITIES_PATH = "/home/utilities"
elif platform.system() == 'Windows':
    UTILITIES_PATH = r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\utilities"
else:
    print('Running on another operating system')


class GapFiller:
    def __init__(self, model_path, universal_model_path, results_path, max_solutions):
        """
        Create a GapFiller object.
        Parameters
        ----------
        model_path
        universal_model_path
        results_path
        """
        self.all_completions = []
        self.optimal_completions = None
        self.all_solutions_with_models = []
        self.positive_solutions = set()
        self.added_demands = None
        self.required_additional_seeds_and_demands = set()
        self.additional_seeds = set()
        self.never_producible = None
        self.reconstructable_targets = None
        self.unproducible = None
        self.resultsPath = results_path
        self.universal_model_copy = None
        self.universal_model_compartmentalized = None
        self.transport_reactions_universal = None
        self.transport_reactions = None
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
        self.reaction_dict = {r.id: r for r in self.universal_model.reactions}
        self.max_solutions = max_solutions

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

    def run(self, optimize: bool = False, write_to_file: bool = False, removed_reactions: list = None, **kwargs):
        """
        Run the gap-filling algorithm.

        Parameters
        ----------
        write_to_file
        removed_reactions
        optimize
        kwargs:
            Keyword arguments to pass to the gap-filling algorithm.

        """
        print("Running gap-filling algorithm...")
        if self.original_model is not None:
            write_sbml_model(self.original_model, self.model_path)
        else:
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
            seeds = sbml.readSBMLseeds(self.seeds_path)
            targets = sbml.readSBMLtargets(self.targets_path)
            print("Getting minimal completion...")
            # Call the get_minimal_completion_size function this order: draftnet, repairnet, seeds, targets)
            from meneco.meneco import query
            model = query.get_unproducible(draftnet, targets, seeds)
            unproducible = set([a[0] for pred in model if pred == 'unproducible_target' for a in model[pred]])
            print(f'Unproducible targets ({len(unproducible)}):\n', '\n'.join(unproducible))
            self.unproducible = unproducible

            repairnet = sbml.readSBMLnetwork(self.universal_model_path, 'repairnet')

            combinet = draftnet
            combinet = TermSet(combinet.union(repairnet))

            model = query.get_unproducible(combinet, targets, seeds)
            self.never_producible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
            print(f'Unreconstructable targets ({len(self.never_producible)}):\n', '\n'.join(self.never_producible))
            print("Identifying additional seeds...")
            if self.never_producible:
                self.identify_additional_seeds(combinet, targets)

            seeds = sbml.readSBMLseeds(self.seeds_path)
            model = query.get_unproducible(combinet, targets, seeds)
            self.never_producible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
            print(f'Unreconstructable targets after additional seeds ({len(self.never_producible)}):\n', '\n'.join(self.never_producible))

            reconstructable_targets = set(list(unproducible.difference(self.never_producible)))
            print(f'Reconstructable targets ({len(reconstructable_targets)}):\n', '\n'.join(reconstructable_targets))
            self.reconstructable_targets = reconstructable_targets

            targets = TermSet(Atom('target("' + t + '")') for t in reconstructable_targets)
            self.minimal_completion = extract_xreactions(get_minimal_completion_size(draftnet, repairnet, seeds, targets), False)
            print("Minimal completion finished.")
            print(self.minimal_completion)
            # the following line does not determine all completions immediately, it only returns a generator to do that.
            self.optimal_completions = query.get_optimal_completions(draftnet, repairnet, seeds, targets, len(self.minimal_completion), nmodels=self.max_solutions)
            print("Optimal completions finished.")

            time_end = time.time()
            self.gftime = time_end - time_start

            print("Time taken: ", self.gftime)
            self.minimal_set = set(
                atom[0] for pred in self.minimal_completion if pred == 'xreaction'
                for atom in self.minimal_completion[pred])

        print("Gap-filling algorithm finished.")

        self.add_reactions_to_model()

        sol = self.original_model.slim_optimize()
        if round(sol, 5) > 0:
            self.remove_unnecessary_seeds_and_demands()
        print(self.original_model.summary())
        write_sbml_model(self.original_model, join(self.resultsPath, "gapfilled_model.xml"))
        return self.report(write_to_file=write_to_file, removed_reactions=removed_reactions)

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
    def from_folder(cls, folder_path, results_path, temporary_universal_model: bool = False, objective_function_id: str = None,
                    compartments=None, max_solutions=1000):
        """
        Create a Gapfilling object from a folder containing the model, seeds, targets and universal model.
        Parameters
        ----------
        folder_path
        results_path
        temporary_universal_model
        objective_function_id
        compartments

        Returns
        -------

        """
        shutil.copy(join(UTILITIES_PATH, "universal_model.xml"), folder_path)
        model_file = join(folder_path, "model.xml")
        seeds_file = join(folder_path, "model_seeds.xml")
        targets_file = join(folder_path, "model_targets.xml")
        universal_model_file = join(folder_path, "universal_model.xml")

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

        model_path = join(folder_path, model_file)
        original_model = Model.from_sbml(model_path, objective_function_id)
        original_model.create_trnas_reactions()

        gap_filler = cls(model_path=model_path, results_path=results_path,
                         universal_model_path=os.path.join(folder_path, universal_model_file),
                         max_solutions=max_solutions)
        gap_filler.original_model = original_model
        gap_filler.seeds_path = os.path.join(folder_path, seeds_file)
        gap_filler.targets_path = os.path.join(folder_path, targets_file)
        gap_filler.clone_model(compartments=compartments, folder_path=folder_path)

        # gap_filler.universal_model_path = os.path.join(folder_path, 'universal_model_compartmentalized.xml')   # to remove, just for debugging
        # gap_filler.universal_model = Model(cobra.io.read_sbml_model(os.path.join(folder_path, 'universal_model_compartmentalized.xml')))

        gap_filler.add_transport_reaction_for_dead_ends(exclude_extr=True, compartments=compartments)

        if temporary_universal_model:
            if not objective_function_id:
                raise ValueError(
                    "Objective function ID must be specified for the creation of a temporary universal model.")
            gap_filler.build_temporary_universal_model(folder_path, True)

            if os.path.isfile(os.path.join(folder_path, 'temporary_universal_model.xml')):
                print('Temporary universal model file successfully created.')
                gap_filler.universal_model_path = os.path.join(folder_path, 'temporary_universal_model.xml')
                print('Temporary universal model file path:', gap_filler.universal_model_path)

                gap_filler.temporary_universal_model = Model(cobra.io.read_sbml_model(
                    os.path.join(folder_path, 'temporary_universal_model.xml')))
                gap_filler.universal_model = gap_filler.temporary_universal_model
            else:
                raise FileNotFoundError("Temporary universal model file not found.")
        return gap_filler

    def add_reactions_to_model(self):
        """
        Add reactions to a model based on the obtained completion.
        For each optimal completion, it adds the reactions. If a positive solution is obtained (fba sol) it adds the solution to a list.
        Otherwise, it adds demands to the model and tries to optimize again. If a positive solution is obtained, it adds the solution to a list.
        In the end, if at least one positive solution is obtained, the first one is selected as the final model. Otherwise, the minimal completion
        is added to the model (even without producing biomass).
        Parameters
        ----------
        """
        positive_solutions = set()
        already_in_original_model = {reaction.id for reaction in self.original_model.reactions}
        seeds = read_sbml_model(self.seeds_path)
        mets_in_medium = set()
        for exchange in self.original_model.exchanges:
            if exchange.lower_bound < 0:
                mets_in_medium.add(exchange.reactants[0].id)

        for m in self.optimal_completions:
            completion = set(a[0] for pred in m if pred == 'xreaction' for a in m[pred])
            self.all_completions.append(list(completion))
            temp_model = self.original_model.copy()
            to_add = []
            for reaction in completion:
                save_reaction = self.universal_model.reactions.get_by_id(reaction.replace("R_", "")).copy()
                if {save_reaction.id}.issubset(already_in_original_model):
                    save_reaction.id = save_reaction.id + "_universal"
                to_add.append(save_reaction)
            model_ids = [e.id for e in temp_model.metabolites]
            for seed in seeds.metabolites:
                if seed.id in model_ids and not {seed.id}.issubset(mets_in_medium):
                    temp_model.create_sink(seed.id)
            already_in_original_model = already_in_original_model.union({reaction.id for reaction in to_add})
            temp_model.add_reactions(to_add)
            sol = temp_model.slim_optimize()
            if round(sol, 5) > 0:
                positive_solutions.add((temp_model, tuple([reaction.id for reaction in to_add])))
            else:
                temp_model, sol, demands = self.add_demands_v2(temp_model)
                if round(sol, 5) > 0:
                    to_add += demands
                    positive_solutions.add((temp_model, tuple([reaction.id for reaction in to_add])))
        self.positive_solutions = positive_solutions
        if self.positive_solutions:
            self.original_model = list(self.positive_solutions)[0][0]
            print(self.original_model.slim_optimize())
        else:
            to_add = []
            for reaction in self.minimal_completion:
                save_reaction = self.universal_model.reactions.get_by_id(reaction.replace("R_", "")).copy()
                if {save_reaction.id}.issubset(already_in_original_model):
                    save_reaction.id = save_reaction.id + "_universal"
                to_add.append(save_reaction)
            model_ids = [e.id for e in self.original_model.metabolites]
            for seed in seeds.metabolites:
                if seed.id in model_ids and not {seed.id}.issubset(mets_in_medium):
                    self.original_model.create_sink(seed.id)
            self.original_model.add_reactions(to_add)
        print("Positive solutions:", positive_solutions)

    def add_demands_v2(self, model):
        dead_ends = self.find_deadends(model)
        demands = []
        for dead_end in dead_ends:
            demands.append(model.create_demand(dead_end.id))
        sol = model.slim_optimize()
        return model, sol, demands

    def report(self, write_to_file=False, removed_reactions: list = None):

        if self.filled_model is None:
            self.filled_model = copy.deepcopy(self.original_model)
        sol = self.original_model.optimize()
        df = self.original_model.summary(sol).to_frame()
        summary_as_list = df.loc[round(df["flux"], 6) != 0].to_dict(orient='records')
        summary_parsed = []
        for row in summary_as_list:
            summary_parsed.append(f"{row}")
        unproducible = list(self.unproducible) if self.unproducible else []
        unreconstructable = list(self.never_producible) if self.never_producible else []
        reconstructable = list(self.reconstructable_targets) if self.reconstructable_targets else []
        minimal_completion = list(self.minimal_completion) if self.minimal_completion else []
        essential_seeds_and_demands = [met.id for met in self.required_additional_seeds_and_demands] if self.required_additional_seeds_and_demands else []
        report_dict = {
            'title': 'Gap-filling report',
            'files_used': {
                'model': os.path.basename(self.model_path),
                'seeds': os.path.basename(self.seeds_path),
                'targets': os.path.basename(self.targets_path),
                'universal_model': os.path.basename(self.universal_model_path)
            },
            'execution_time': self.gftime if hasattr(self, 'gftime') else None,
            'Objective': sol.objective_value,
            'artificial_removed_reactions': removed_reactions if removed_reactions else [],
            'unproducible_targets': unproducible,
            'unreconstructable_targets': unreconstructable,
            'reconstructable_targets': reconstructable,
            'minimal_completion': minimal_completion,
            'additional_seeds': [met.id for met in self.additional_seeds] if self.additional_seeds else [],
            'essential_additional_seeds_and_demands': essential_seeds_and_demands,
            'summary': summary_parsed,
            'positive_solutions': [solution[1] for solution in self.positive_solutions] if self.positive_solutions else [],
            'all_completions': self.all_completions
        }
        print(report_dict)
        if write_to_file:
            json_report_path = os.path.join(self.resultsPath, "gapfilling_report.json")
            with open(json_report_path, 'w') as json_file:
                json.dump(report_dict, json_file, indent=4)

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

        self.universal_model_copy = model_to_use.copy()

        self.universal_model_compartmentalized = cobra.Model()

        self.clone_metabolites(compartments)

        self.clone_reactions(compartments)

        self.clone_groups(compartments)

        self.universal_model = self.universal_model_compartmentalized
        self.universal_model_path = os.path.join(kwargs['folder_path'], 'universal_model_compartmentalized.xml')

        # get the transport reactions
        # cloned_transport_reactions = self.get_transport_reactions_of_draft(False)
        # for reaction in self.universal_model_compartmentalized.reactions:
        #         reaction.bounds = (-1000, 1000)
        demands = []
        for metabolite in self.universal_model_compartmentalized.metabolites:
            if len(metabolite.reactions) <= 1:
                # add demand reaction
                demand_reaction = cobra.Reaction('DM_' + metabolite.id, name="Demand reaction for " + metabolite.id)
                demand_reaction.add_metabolites({metabolite: -1})
                demand_reaction.bounds = (0, 1000)
                demands.append(demand_reaction)
        self.universal_model_compartmentalized.add_reactions(demands)
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
                    metabolite_in_compartment = self.universal_model_copy.metabolites.get_by_id(
                        new_metabolite_id_in_compartment)
                    cloned_reaction.add_metabolites({metabolite: -st, metabolite_in_compartment: st})
                cloned_reactions.append(cloned_reaction)
        return cloned_reactions

    def clone_metabolites(self, compartments):
        parallelbar.parallelbar.PROGRESS = "CLONING METABOLITES"
        list_of_batches = [self.universal_model_copy.metabolites[i:i + 100] for i in
                           range(0, len(self.universal_model_copy.metabolites), 100)]
        new_metabolites = progress_imap(partial(clone_metabolite, compartments), list_of_batches)
        new_metabolites = [item for sublist in new_metabolites for item in sublist]
        self.universal_model_copy.add_metabolites(new_metabolites)

    def get_transport_reactions_of_draft(self, clone: bool = False) -> list:
        """
        Get all the transport reactions in the draft model
        """
        self.transport_reactions = [r for r in self.original_model.reactions if len(r.compartments) > 1]

        if clone:
            return [rxn.copy() for rxn in self.transport_reactions]
        else:
            return self.transport_reactions

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
        solver
        objective
        objective_function_id
        """

        if self.objective_function_id is None:
            self.objective_function_id = objective_function_id

        # Set the solver
        set_solver(self, solver)

        # Run BioISO
        bio = BioISO(self.objective_function_id, self.original_model, objective)
        bio.run(2, False)

        # Get the results
        results = bio.get_tree()
        biomass_components = results["M_root_M_root_M_root_product"]["next"]

        dead_ends = []

        for biomass_component in biomass_components:
            biomass_component_role = biomass_components[biomass_component].get("role")
            specific_biomass_components = biomass_components[biomass_component].get("next")
            for specific_biomass_component in specific_biomass_components:
                specific_biomass_component_report = specific_biomass_components[specific_biomass_component]
                analysis = specific_biomass_component_report.get("analysis")
                role = specific_biomass_component_report.get("role")

                if not analysis and ((role == "Reactant" and biomass_component_role == "Reactant") or
                                     (role == "Product" and biomass_component_role == "Product")):
                    metabolite_id = specific_biomass_component_report.get("identifier")
                    if metabolite_id not in dead_ends:
                        dead_ends.append(metabolite_id)

        self.dead_ends = dead_ends

        return dead_ends

    def add_known_transport_reactions(self):
        """
        Add known transport reactions to the model to resolve dead-end metabolites.
        """
        if self.dead_ends is None:
            dead_ends = self.get_dead_end_metabolites()
        else:
            dead_ends = self.dead_ends

        transport_reactions = {met: [] for met in dead_ends}
        for reaction in self.transport_reactions:
            for met in reaction.metabolites:
                if met in transport_reactions:
                    transport_reactions[met].append(reaction.id)

        for met in dead_ends:
            for reaction_id in transport_reactions[met]:
                reaction = self.reaction_dict[reaction_id]
                self.original_model.add_reactions([reaction.copy()])

        print(f"Added known transport reactions for dead-end metabolites.")

    def clone_reactions(self, compartments):
        """
        Clone reactions in the model.
        Parameters
        ----------
        compartments

        Returns
        -------

        """
        parallelbar.parallelbar.PROGRESS = "CLONING REACTIONS"
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

    @staticmethod
    def met_is_being_consumed_or_produced_at_reaction(met_id: str, reaction: cobra.Reaction):
        """
        This function checks if a metabolite is being consumed or produced at a reaction

        Parameters
        ----------
        met_id: str
            The ID of the metabolite.
        reaction: cobra.Reaction
            The reaction to check.

        Returns
        -------
        str
            'consumed' if the metabolite is a reactant in the reaction,
            'produced' if the metabolite is a product of the reaction,
            'not involved' if the metabolite is neither a reactant nor a product of the reaction.
        """
        metabolite = reaction.model.metabolites.get_by_id(met_id)
        if metabolite in reaction.reactants:
            return 'consumed'
        elif metabolite in reaction.products:
            return 'produced'
        else:
            return 'not involved'

        # for met in dead_end_metabolites:...

    def check_metabolite_in_reaction(self, met_id: str, reactions: list = None) -> dict:
        """
        This function checks if a metabolite is being consumed or produced at a reaction

        Parameters
        ----------
        met_id: str
            The ID of the metabolite.
        reactions: list[cobra.Reaction]
            The reaction to check.

        Returns
        -------
        str
            'consumed' if the metabolite is a reactant in the reaction,
            'produced' if the metabolite is a product of the reaction,
            'not involved' if the metabolite is neither a reactant nor a product of the reaction.
        """
        if reactions is None:
            reactions = self.original_model.metabolites.get_by_id(met_id).reactions

        involvement = {}
        i = 0
        for reaction in reactions:
            involvement[reaction.id] = self.met_is_being_consumed_or_produced_at_reaction(met_id, reaction)
            print('reaction', i, ': ', reaction.id,
                  'is reversible: ', reaction.reversibility,
                  'and the lower bound is: ', reaction.lower_bound,
                  'and the upper bound is: ', reaction.upper_bound,
                  'and the objective coefficient is: ', reaction.objective_coefficient)
            i += 1

        return involvement

    def search_dead_end(self, metabolite_id, consumed=True, produced=True, verbose=False) -> tuple:
        """
        This method searches for reactions in the universal model where a given metabolite is consumed or produced.

        Parameters
        ----------
        verbose: bool, optional
            If True, the method will print the reactions where the metabolite is consumed or produced.
        metabolite_id: str
            The ID of the metabolite.
        consumed: bool, optional
            If True, the method will search for reactions where the metabolite is consumed.
        produced: bool, optional
            If True, the method will search for reactions where the metabolite is produced.

        Returns
        -------
        tuple
            A tuple containing two lists:
            - The first list contains reactions where the metabolite is consumed.
            - The second list contains reactions where the metabolite is produced.
        """
        # Get the metabolite object from the model
        metabolite = self.universal_model_compartmentalized.metabolites.get_by_id(metabolite_id)

        # Initialize lists to store reactions where the metabolite is consumed or produced
        consumed_reactions = []
        produced_reactions = []

        # Iterate over all reactions in the model
        for reaction in self.universal_model_compartmentalized.reactions:
            # Check if the metabolite is a reactant in the reaction
            if metabolite in reaction.reactants and consumed:
                consumed_reactions.append(reaction)
            # Check if the metabolite is a product in the reaction
            if metabolite in reaction.products and produced:
                produced_reactions.append(reaction)

        if verbose:
            if consumed:
                print(f"Reactions where {metabolite_id} is consumed:")
                for reaction in consumed_reactions:
                    print(reaction.id)
            if produced:
                print(f"Reactions where {metabolite_id} is produced:")
                for reaction in produced_reactions:
                    print(reaction.id)

        return consumed_reactions, produced_reactions

    def deal_dead_ends(self, dead_ends: list):
        """
        This method deals with dead-end metabolites by searching for reactions in the universal model where the dead-end metabolite is consumed or produced and then adding those reactions to the draft model.

        Parameters
        ----------
        dead_ends: list
            List of dead-end metabolite IDs.

        Returns
        -------
        cobra.Model
            The updated draft model with added reactions.
        """

        if dead_ends is None:
            dead_ends = self.dead_ends
            if dead_ends is None:
                self.identify_dead_end_metabolites_with_bioiso(objective_function_id=self.objective_function_id)
                dead_ends = self.dead_ends

        # Copy the draft model
        draft = self.original_model.copy()

        # Iterate over dead ends
        for dead in dead_ends:
            # Reset flags
            flag_consumed = False
            flag_produced = False

            # Initialize lists to store reactions where the dead end is consumed or produced
            consumed_uni = []
            produced_uni = []

            # Find reactions associated with dead ends in the draft model
            reactions_dead_end = draft.metabolites.get_by_id(dead).reactions

            # Check if dead end is being produced or consumed
            for reaction in reactions_dead_end:
                if self.met_is_being_consumed_or_produced_at_reaction(dead, reaction) == 'consumed':
                    flag_consumed = True
                elif self.met_is_being_consumed_or_produced_at_reaction(dead, reaction) == 'produced':
                    flag_produced = True

            # Find reactions in the universal model where the dead end is consumed or produced
            if flag_consumed and flag_produced:
                consumed_uni, produced_uni = self.search_dead_end(dead, consumed=True, produced=True)
            elif flag_consumed:
                consumed_uni, produced_uni = self.search_dead_end(dead, consumed=False, produced=True)
            elif flag_produced:
                consumed_uni, produced_uni = self.search_dead_end(dead, consumed=True, produced=False)

            # Add reactions from the universal model to the draft model
            for reaction in consumed_uni + produced_uni:
                # Create a copy of the reaction
                reaction_copy = reaction.copy()

                # Check if the metabolite is in different compartments
                if dead.split('__')[-1] not in [met.id.split('__')[-1] for met in reaction_copy.metabolites]:
                    # Create a new transport reaction
                    transport_reaction = Reaction(id=f"transport_{dead}")
                    transport_reaction.name = f"Transport of {dead}"
                    transport_reaction.lower_bound = -1000
                    transport_reaction.upper_bound = 1000

                    # Create two metabolite objects for each compartment # review here
                    met1 = Metabolite(id=f"{dead}", compartment=list(reaction_copy.compartments)[0])
                    met2 = Metabolite(id=f"{dead}", compartment=list(reaction_copy.compartments)[1])

                    # Add the metabolites to the transport reaction
                    transport_reaction.add_metabolites({met1: -1, met2: 1})

                    # Add the transport reaction to the draft model
                    draft.add_reactions([transport_reaction])

                # Add the reaction to the draft model
                draft.add_reactions([reaction_copy])

        return draft

    def add_transport_reaction_for_dead_ends(self, compartments: list = None, exclude_extr: bool = False):
        parallelbar.parallelbar.PROGRESS = "ADDING TRANSPORT"
        dead_ends = cobra.io.read_sbml_model(self.targets_path).metabolites

        # dead_ends = self.universal_model.metabolites

        if compartments is None:
            compartments = get_compartments_in_model(self.original_model)

        if exclude_extr:
            for extr_id in ['extr', 'e', 'extracellular']:
                if extr_id in compartments:
                    compartments.remove(extr_id)

        if len(compartments) < 2:
            print("Warning: There must be at least two compartments to create a transport reaction.")
            return None

        compartments_combinations = list(itertools.combinations(compartments, 1))

        compartments_combinations = list(set([(comp1[0], 'cytop') for comp1 in compartments_combinations if comp1[0] != 'cytop']))

        metabolites_ids = list(set(['__'.join(met.id.split("__")[:-1]) for met in dead_ends]))

        list_of_batches = [metabolites_ids[i:i + 500] for i in range(0, len(metabolites_ids), 500)]

        result = progress_imap(partial(self.create_transport, compartments_combinations), list_of_batches)
        for batch in result:
            self.universal_model_compartmentalized.add_reactions(batch)

    def create_transport(self, compartments_combinations, metabolites):
        to_add = []
        for metabolite_id in metabolites:
            for compartments_combination in compartments_combinations:
                transport_reaction = Reaction(id=f"T_{metabolite_id}_{compartments_combination[0]}_to_{compartments_combination[1]}")
                transport_reaction.name = f"Transport of {metabolite_id} between {compartments_combination[0]} and {compartments_combination[1]}"
                transport_reaction.lower_bound = -1000
                transport_reaction.upper_bound = 1000
                met1 = self.universal_model.metabolites.get_by_id(f"{metabolite_id}__{compartments_combination[0]}")
                met2 = self.universal_model.metabolites.get_by_id(f"{metabolite_id}__{compartments_combination[1]}")
                transport_reaction.add_metabolites({met1: -1, met2: 1})
                to_add.append(transport_reaction)
        return to_add

    def build_temporary_universal_model(self, folder_path, related_pathways: bool = False):
        """
        Build a temporary universal model from the universal model and the gap-filling results.
        Parameters
        ----------
        related_pathways
        folder_path

        Returns
        -------

        """
        pathways_to_ignore = {'Metabolic pathways', 'Biosynthesis of secondary metabolites', 'Microbial metabolism in diverse environments',
                              'Biosynthesis of cofactors', 'Carbon metabolism', 'Fatty acid metabolism'}
        pathways_to_keep = ["Nicotinate and nicotinamide metabolism", "One carbon pool by folate, "
                              "Folate biosynthesis","Pantothenate and CoA biosynthesis", "Riboflavin metabolism"]
        metabolite_pathways_map = {}
        universal_model = Model(model=self.universal_model_compartmentalized)
        print('Number of reactions in universal model:', len(universal_model.reactions))
        targets = read_sbml_model(self.targets_path)
        cofactors = json.load(open(os.path.join(UTILITIES_PATH, 'cofactors.json'), 'r'))
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
            if os.path.exists(os.path.join(UTILITIES_PATH, 'related_pathways_map.json')):
                with open(os.path.join(UTILITIES_PATH, 'related_pathways_map.json'), 'r') as related_pathways_map_file:
                    related_pathways_map = json.load(related_pathways_map_file)
            else:
                related_pathways_map = create_related_pathways_map(universal_model, folder_path)

            for pathway in pathways_to_keep:
                related_pathways.update(get_related_pathways(pathway, related_pathways_map))
            pathways_to_keep = pathways_to_keep.union(related_pathways) - pathways_to_ignore

        for metabolite, pathways in metabolite_pathways_map.items():
            pathways_copy = pathways.copy()
            for pathway in pathways_copy:
                metabolite_pathways_map[metabolite].update(get_related_pathways(pathway, related_pathways_map))
        to_keep = set([reaction for reaction in universal_model.reactions if "T_" in reaction.id])
        for pathway in universal_model.groups:
            if pathway.name in pathways_to_keep:
                to_keep.update(reaction for reaction in pathway.members)
        for reaction_id in universal_model.reaction_pathway_map.keys():
            if len(universal_model.reaction_pathway_map[reaction_id]) == 0:
                to_keep.add(universal_model.reactions.get_by_id(reaction_id))
        groups = [pathway for pathway in universal_model.groups if pathway.name in pathways_to_keep]
        self.universal_model_compartmentalized = Model()
        self.universal_model_compartmentalized.add_reactions(to_keep)
        self.universal_model_compartmentalized.add_groups(groups)
        print('Number of reactions in temporary universal model:', len(self.universal_model_compartmentalized.reactions))
        write_sbml_model(self.universal_model_compartmentalized, os.path.join(folder_path, 'temporary_universal_model.xml'))

    def identify_additional_seeds(self, combinet, targets):
        cofactors = json.load(open(os.path.join(UTILITIES_PATH, 'cofactors.json'), 'r'))
        for metabolite in self.never_producible:
            metabolite_id = metabolite.lstrip("M_")
            if metabolite_id in self.universal_model.metabolite_pathway_map.keys():
                for pathway in self.universal_model.metabolite_pathway_map[metabolite_id]:
                    if pathway in self.universal_model.pathway_metabolites_map.keys():
                        metabolites_in_pathway = self.universal_model.pathway_metabolites_map[pathway]
                        for potential_cofactor in metabolites_in_pathway:
                            potential_cofactor_id = '__'.join(potential_cofactor.split("__")[:-1])
                            if potential_cofactor_id in cofactors.keys():
                                self.additional_seeds.add(self.universal_model.metabolites.get_by_id(potential_cofactor))
        current_seeds = set(read_sbml_model(self.seeds_path).metabolites)
        final_seeds = set(current_seeds.union(self.additional_seeds))
        res = []
        for met in final_seeds:
            compartment = '__'.join(met.id.split("__")[-1:])
            res.append((met.id, compartment))

        write_metabolites_to_sbml("model_seeds.xml", dirname(self.seeds_path), list(res))

        # first, I tried to remove the additional seeds that are not necessary. However, even if a seed is not necessary, not having them increases the running time a lot.

        # seeds = sbml.readSBMLseeds(self.seeds_path)
        # model = query.get_unproducible(combinet, targets, seeds)
        # never_producible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
        # to_remove = set()
        # for seed in sorted(self.additional_seeds, key=lambda x: x.id, reverse=True):
        #     old_number_of_never_producible = len(never_producible)
        #     res.pop(res.index((seed.id, '__'.join(seed.id.split("__")[-1:]))))
        #     write_metabolites_to_sbml("model_seeds.xml", dirname(self.seeds_path), list(res))
        #     seeds = sbml.readSBMLseeds(self.seeds_path)
        #     model = query.get_unproducible(combinet, targets, seeds)
        #     never_producible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
        #     if len(never_producible)>old_number_of_never_producible:
        #          res.append((seed.id, '__'.join(seed.id.split("__")[-1:])))
        #     else:
        #         to_remove.add(seed)
        # self.additional_seeds = self.additional_seeds - to_remove
        print(f"Number of additional seeds: {len(self.additional_seeds)}")
        print(f"Additional seeds: {[met.id for met in self.additional_seeds]}")
        write_metabolites_to_sbml("model_seeds.xml", dirname(self.seeds_path), list(res))

    @deprecated
    def add_demands(self):
        for solution, temp_model in self.all_solutions_with_models:
            dead_ends = self.find_deadends(temp_model)
            to_remove = set()
            for dead_end in dead_ends:
                temp_model.create_demand(dead_end.id)
            demands = {reaction for reaction in self.original_model.reactions if reaction.id.startswith("DM_")}
            initial_sol = self.original_model.slim_optimize()
            if round(initial_sol, 5) > 0:
                for demand in demands:
                    with temp_model as m:
                        m.reactions.get_by_id(demand.id).bounds = (0, 0)
                        sol = m.slim_optimize()
                        if round(sol, 5) > 0:
                            to_remove.add(demand)
                        else:
                            m.reactions.get_by_id(demand.id).bounds = (0, 1000)
            self.original_model.remove_reactions(list(to_remove))
            self.added_demands = demands - to_remove
            sol = self.original_model.slim_optimize()
            if round(sol, 5) > 0:
                solution.add(self.added_demands)
                self.positive_solutions.add(solution)
        if self.positive_solutions:
            self.original_model.add_reactions(self.universal_model.reactions.get_by_id(rxn.id) for rxn in list(self.positive_solutions)[0])

    def find_deadends(self, model=None):
        """
        Return metabolites that are only produced in reactions.

        Metabolites that are involved in an exchange reaction are never
        considered to be dead ends.

        """

        if model is None:
            model = self.original_model

        def is_only_product(metabolite: cobra.Metabolite, reaction: cobra.Reaction) -> bool:
            """Determine if a metabolite is only a product of a reaction."""
            if reaction.reversibility:
                return False
            if reaction.get_coefficient(metabolite) > 0:
                return reaction.lower_bound >= 0.0 and reaction.upper_bound > 0
            else:
                return reaction.lower_bound < 0.0 and reaction.upper_bound <= 0

        exchanges = frozenset(model.exchanges)
        return [
            met
            for met in model.metabolites
            if (len(met.reactions) > 0)
               and all((rxn not in exchanges) and is_only_product(met, rxn) for rxn in met.reactions)
        ]

    def remove_unnecessary_seeds_and_demands(self):
        """
        Remove unnecessary seeds and demands from the model.
        Returns
        -------

        """
        self.required_additional_seeds_and_demands = self.additional_seeds.union(set([met.id for demand in self.original_model.demands for met in demand.metabolites]))
        to_remove = set()
        for reaction in self.original_model.reactions:
            if reaction.id.startswith("Sk_") or reaction.id.startswith("DM_"):
                original_bounds = reaction.bounds
                reaction.bounds = (0, 0)
                sol = self.original_model.slim_optimize()
                if round(sol, 3) > 0:
                    to_remove.add(reaction)
                else:
                    reaction.bounds = original_bounds
        self.original_model.remove_reactions(list(to_remove))
        print(to_remove)
        for sink in to_remove:
            self.required_additional_seeds_and_demands = self.required_additional_seeds_and_demands - set([met.id for met in sink.metabolites])
