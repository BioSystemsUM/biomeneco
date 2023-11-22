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
from .utils import write_metabolites_to_sbml
from gap_filling_dl.kegg_api import create_related_pathways_map, get_related_pathways
from .model import Model
from .utils import get_compartments_in_model, clone_metabolite
from .utils import identify_dead_end_metabolites

### The following code is to set up the progress bar. The parallelbar package does not allow changing the description is
# the bar, so I had to change the code.
# Set the default progress bar description to "DONE"
parallelbar.parallelbar.PROGRESS = "DONE"


def _my_process_status(bar_size, worker_queue):
    """
        Custom process status function to update the progress bar and error bar based on the worker queue.

        This function overrides the default behavior to allow for a secondary 'error bar' which
        tracks the number of errors encountered during processing. The main progress bar tracks
        successful updates.

        Args:
            bar_size (int): The total size of the progress bar, representing the work to be done.
            worker_queue (queue.Queue): A queue that receives updates from worker threads/processes.

        The function listens for updates from the worker queue and updates the progress bar accordingly.
        If an error is encountered (indicated by a True flag in the queue), an error bar is updated.
        When a None flag is received, it signals the end of processing, and the function closes the
        progress bars and exits.
        """
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
elif platform.system() == 'Darwin':
    UTILITIES_PATH = "/Users/josediogomoura/gap_filling_dl/utilities"
else:
    print('Running on another operating system')


class GapFiller:
    def __init__(self, model_path: str, universal_model_path: str, results_path: str, max_solutions: int):
        """
        Initialize the GapFiller object which is designed to fill gaps in metabolic models.

        The GapFiller uses a given metabolic model and a universal model containing all known reactions to find solutions that enable the production of target metabolites. It can generate multiple solutions and identify the minimal set of reactions needed to achieve the desired metabolic functions.

        Parameters
        ----------
        model_path : str
            Path to the file containing the metabolic model to be gap-filled.
        universal_model_path : str
            Path to the file containing the universal model with all known reactions.
        results_path : str
            Path where the results of the gap-filling process will be stored.
        max_solutions : int
            Maximum number of solutions to find during the gap-filling process.

        Attributes
        ----------
        all_completions : list
            to store all possible completions found during the gap-filling process.
        optimal_completions : list or None
            List to store the optimal completions if identified, otherwise None.
        all_solutions_with_models : list
            to store all solutions along with their corresponding models.
        positive_solutions : set
            to store solutions that result in a positive outcome, such as growth or production of target metabolites.
        added_demands : set or None
            Set to keep track of demand reactions added to the model, or None if not applicable.
        required_additional_seeds_and_demands : set
            to keep track of additional seeds and demands required for the model to produce target metabolites.
        additional_seeds : set
            of additional seed metabolites identified during the gap-filling process.
        never_producible : set or None
            Set of metabolites that are identified as never producible, or None if not applicable.
        reconstructable_targets : set or None
            Set of targets that can be reconstructed, or None if not applicable.
        unproducible : set or None
            Set of metabolites that cannot be produced, or None if not applicable.
        resultsPath : str
            Path to the results directory.
        universal_model_copy : cobra.Model or None
            A copy of the universal model, or None if not yet created.
        universal_model_compartmentalized : cobra.Model or None
            The universal model with compartmentalization, or None if not yet created.
        transport_reactions_universal : list or None
            List of transport reactions in the universal model, or None if not applicable.
        transport_reactions : list or None
            List of transport reactions added to the model, or None if not applicable.
        dead_ends : list or None
            List of dead-end metabolites in the model, or None if not yet identified.
        cloned_model : cobra.Model or None
            A cloned version of the original model, or None if not yet created.
        temporary_universal_model : cobra.Model or None
            A temporary universal model used during the gap-filling process, or None if not yet created.
        original_model : cobra.Model or None
            The original metabolic model before gap-filling, or None if not yet loaded.
        _universal_model : cobra.Model or None
            Internal attribute for the universal model, or None if not yet loaded.
        minimal_set : list or None
            List of reactions constituting the minimal set required for gap-filling, or None if not yet identified.
        minimal_completion : list or None
            List of reactions constituting the minimal completion, or None if not yet identified.
        gftime : float or None
            Time taken for the gap-filling process, or None if not yet computed.
        cobra_filled_model : cobra.Model or None
            The gap-filled model using COBRApy, or None if not yet created.
        results_meneco : dict or None
            Results from the Meneco tool, or None if not yet computed.
        universal_model_path : str
            Path to the universal model file.
        model_path : str
            Path to the metabolic model file.
        _seeds_path : str or None
            Internal path to the seeds file, or None if not yet set.
        _targets_path : str or None
            Internal path to the targets file, or None if not yet set.
        filled_model : cobra.Model or None
            The final gap-filled model, or None if not yet created.
        objective_function_id : str or None
            ID of the objective function used in the model, or None if not yet set.
        reaction_dict : dict
            Dictionary mapping reaction IDs to reaction objects in the universal model.
        max_solutions : int
            The maximum number of gap-filling solutions to compute.
        """

        # Initialize all attributes
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
        """
        Get the file path to the seeds used for gap-filling.

        Returns:
            str: The file path to the seeds.
        """
        return self._seeds_path

    @seeds_path.setter
    def seeds_path(self, file_path):
        """
        Set the file path to the seeds used for gap-filling.

        Args:
            file_path (str): The file path to the seeds.
        """
        self._seeds_path = file_path

    @property
    def targets_path(self):
        """
        Get the file path to the targets used for gap-filling.

        Returns:
            str: The file path to the targets.
        """
        return self._targets_path

    @targets_path.setter
    def targets_path(self, file_path):
        """
        Set the file path to the targets used for gap-filling.

        Args:
            file_path (str): The file path to the targets.
        """
        self._targets_path = file_path

    # read the universal model: cobra.io.read_sbml_model(universal_model_path)

    @property
    def universal_model_path(self):
        """
        Get the file path to the universal model used for gap-filling.

        If `_universal_model_path` is `None`, it sets a default path based on the `model_path`.

        Returns:
            str: The file path to the universal model.
        """
        if self._universal_model_path is None:
            self._universal_model_path = os.path.join(os.path.dirname(self.model_path), 'universal_model.xml')
        return self._universal_model_path

    @universal_model_path.setter
    def universal_model_path(self, value):
        """
        Set the file path to the universal model used for gap-filling.

        Args:
            value (str): The file path to the universal model.
        """
        self._universal_model_path = value

    @property
    def universal_model(self):
        """
        Get the universal model used for gap-filling.

        If `_universal_model` is `None`, it reads the universal model from the `universal_model_path`.

        Returns:
            cobra.Model: The universal model.
        """
        if self._universal_model is None:
            self._universal_model = cobra.io.read_sbml_model(self.universal_model_path)
        return self._universal_model

    @universal_model.setter
    def universal_model(self, value):
        """
            Set the universal model used for gap-filling.

            Args:
                value (cobra.Model): The universal model to set.
            """
        self._universal_model = value

    @property
    def original_model(self):
        """
            Get the original model used for gap-filling.

            If `_original_model` is `None`, it reads the original model from the `model_path`.

            Returns:
                cobra.Model: The original model.
            """
        if self._original_model is None:
            self._original_model = cobra.io.read_sbml_model(self.model_path)
        return self._original_model

    @original_model.setter
    def original_model(self, value):
        """
            Set the original model used for gap-filling.

            Args:
                value (cobra.Model): The original model to set.
            """
        self._original_model = value

    def run(self, optimize: bool = False, write_to_file: bool = False, removed_reactions: list = None, **kwargs):
        """
        Run the gap-filling algorithm.

        Parameters
        ----------
        optimize : bool, optional
            If True, the gap-filling algorithm will optimize the model after adding reactions.
        write_to_file : bool, optional
            If True, the gap-filled model will be written to a file.
        removed_reactions : list, optional
            List of reactions to remove from the model before gap-filling.
        kwargs
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
            print(f'Unreconstructable targets after additional seeds ({len(self.never_producible)}):\n',
                  '\n'.join(self.never_producible))

            reconstructable_targets = set(list(unproducible.difference(self.never_producible)))
            print(f'Reconstructable targets ({len(reconstructable_targets)}):\n', '\n'.join(reconstructable_targets))
            self.reconstructable_targets = reconstructable_targets

            targets = TermSet(Atom('target("' + t + '")') for t in reconstructable_targets)
            self.minimal_completion = extract_xreactions(
                get_minimal_completion_size(draftnet, repairnet, seeds, targets), False)
            print("Minimal completion finished.")
            print(self.minimal_completion)
            # the following line does not determine all completions immediately, it only returns a generator to do that.
            self.optimal_completions = query.get_optimal_completions(draftnet, repairnet, seeds, targets,
                                                                     len(self.minimal_completion),
                                                                     nmodels=self.max_solutions)
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
    def from_folder(cls, folder_path, results_path, temporary_universal_model: bool = False,
                    objective_function_id: str = None,
                    compartments=None, max_solutions=1000):
        """
        Create a Gapfilling object from a folder containing the model, seeds, targets and universal model.
        Parameters
        ----------
        folder_path : str
            The path to the folder containing the model, seeds, targets and universal model.
        results_path : str
            The path to the folder where the results will be stored.
        temporary_universal_model : bool, optional
            If True, a temporary universal model will be created for the gap-filling process.
        objective_function_id : str, optional
            The ID of the objective function to use for the gap-filling process.
        compartments : list, optional
            A list of compartments to use for the gap-filling process.
        max_solutions : int, optional
            The maximum number of solutions to compute.
        Returns
        -------
        gap_filler : GapFiller
            The GapFiller object.

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
        positive_solutions = set()  # those that produce biomass
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
                temp_model, sol, demands = self.add_demands_v2(temp_model)  # type: ignore
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
                    self.original_model.create_sink(seed.id)  # type: ignore
            self.original_model.add_reactions(to_add)
        print("Positive solutions:", positive_solutions)

    def add_demands_v2(self, model):
        """
        Add demands to the model.

        Parameters
        ----------
        model : cobra.Model
            The model to add demands to.

        Returns
        -------
        model : cobra.Model
            The model with demands added.
        sol : float
            The solution of the model.
        demands : list
            The list of demands added to the model.
        """

        dead_ends = self.find_deadends(model)
        demands = []
        for dead_end in dead_ends:
            demands.append(model.create_demand(dead_end.id))
        sol = model.slim_optimize()
        return model, sol, demands

    def report(self, write_to_file=False, removed_reactions: list = None):
        """
        Generate a report of the gap-filling process.
        Parameters
        ----------
        write_to_file : bool, optional
            If True, the report will be written to a file.
        removed_reactions : list, optional
            List of reactions removed from the model before gap-filling.

        Returns
        -------
        report_dict : dict
            A dictionary containing information about the gap-filling process.

        """

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
        essential_seeds_and_demands = [met.id for met in
                                       self.required_additional_seeds_and_demands] if self.required_additional_seeds_and_demands else []
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
            'positive_solutions': [solution[1] for solution in
                                   self.positive_solutions] if self.positive_solutions else [],
            'all_completions': self.all_completions
        }
        print(report_dict)
        if write_to_file:
            json_report_path = os.path.join(self.resultsPath, "gapfilling_report.json")
            with open(json_report_path, 'w') as json_file:
                json.dump(report_dict, json_file, indent=4)

    def clone_model(self, compartments: list = None, **kwargs):
        """
        Clone the original model and add the compartments to the metabolites and reactions.
        Parameters
        ----------
        compartments : list, optional
            List of compartments to add to the model.

        kwargs :
            Keyword arguments to pass to the compartmentalization function.

        Returns
        -------

        """
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
        """
        Clone a reaction and add the compartments to the metabolites.
        Parameters
        ----------
        compartments : list
            List of compartments to add to the reaction.
        reactions : list
            List of reactions to clone.

        Returns
        -------

        """
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
        """
        Clone metabolites and add the compartments to the metabolites.
        Parameters
        ----------
        compartments : list
            List of compartments to add to the metabolites.

        Returns
        -------
        """
        parallelbar.parallelbar.PROGRESS = "CLONING METABOLITES"
        list_of_batches = [self.universal_model_copy.metabolites[i:i + 100] for i in
                           range(0, len(self.universal_model_copy.metabolites), 100)]
        new_metabolites = progress_imap(partial(clone_metabolite, compartments), list_of_batches)
        new_metabolites = [item for sublist in new_metabolites for item in sublist]
        self.universal_model_copy.add_metabolites(new_metabolites)

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
        solver : str, optional
            The solver to use. Either 'cplex' or 'glpk'.
        objective : str, optional
            The type of objective function. Either 'maximize' or 'minimize'.
        objective_function_id : str
            The ID of the objective function to use for the gap-filling process.

        Returns
        -------
        dead_ends : list
            List of dead-end metabolites.
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

    def clone_reactions(self, compartments):
        """
        Clone reactions in the model.
        Parameters
        ----------
        compartments : list
            List of compartments to add to the reactions.

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
        Parameters
        ----------
        compartments : list
            List of compartments to add to the groups.
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

    def add_transport_reaction_for_dead_ends(self, compartments: list = None, exclude_extr: bool = False,
                                             ):
        """
        Add transport reactions for dead-end metabolites.
        Parameters
        ----------
        mock_model : cobra.Model, optional
            If True, a mock model will be used to identify dead-end metabolites. ONLY USE THIS FOR TESTING.
        compartments : list, optional
            List of compartments to add to the transport reactions.
        exclude_extr : bool, optional
            If True, exclude the extracellular compartment.

        Returns
        -------
        """

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

        compartments_combinations = list(
            set([(comp1[0], 'cytop') for comp1 in compartments_combinations if comp1[0] != 'cytop']))

        metabolites_ids = list(set(['__'.join(met.id.split("__")[:-1]) for met in dead_ends]))

        list_of_batches = [metabolites_ids[i:i + 500] for i in range(0, len(metabolites_ids), 500)]

        result = progress_imap(partial(self.create_transport, compartments_combinations), list_of_batches)
        for batch in result:
            self.universal_model_compartmentalized.add_reactions(batch)

    def create_transport(self, compartments_combinations, metabolites):
        """
        Create transport reactions for metabolites.
        Parameters
        ----------
        compartments_combinations : list
            List of compartments combinations.
        metabolites : list
            List of metabolites to create transport reactions for.

        Returns
        -------
        to_add : list
            List of transport reactions to add to the model.

        """
        to_add = []
        for metabolite_id in metabolites:
            for compartments_combination in compartments_combinations:
                transport_reaction = Reaction(
                    id=f"T_{metabolite_id}_{compartments_combination[0]}_to_{compartments_combination[1]}")
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
        related_pathways : bool, optional
            If True, related pathways will be added to the model.
        folder_path : str
            Path to the folder where the temporary universal model will be saved.

        Returns
        -------

        """
        pathways_to_ignore = {'Metabolic pathways', 'Biosynthesis of secondary metabolites',
                              'Microbial metabolism in diverse environments',
                              'Biosynthesis of cofactors', 'Carbon metabolism', 'Fatty acid metabolism'}
        pathways_to_keep = ["Nicotinate and nicotinamide metabolism", "One carbon pool by folate, "
                                                                      "Folate biosynthesis",
                            "Pantothenate and CoA biosynthesis", "Riboflavin metabolism"]
        metabolite_pathways_map = {}
        universal_model = Model(model=self.universal_model_compartmentalized)
        print('Number of reactions in universal model:', len(universal_model.reactions))
        targets = read_sbml_model(self.targets_path)
        cofactors = json.load(open(os.path.join(UTILITIES_PATH, 'cofactors.json'), 'r'))
        for target in targets.metabolites:
            if target.id in universal_model.metabolite_pathway_map.keys():
                if '__'.join(target.id.split("__")[:-1]) in cofactors.keys():
                    pathways_to_keep += [pathway for pathway in universal_model.metabolite_pathway_map[target.id] if
                                         pathway in cofactors['__'.join(target.id.split("__")[:-1])]]
                    metabolite_pathways_map[target.id] = set(
                        pathway for pathway in universal_model.metabolite_pathway_map[target.id] if
                        pathway in cofactors['__'.join(target.id.split("__")[:-1])])
                else:
                    pathways_to_keep += [pathway for pathway in universal_model.metabolite_pathway_map[target.id]]
                    metabolite_pathways_map[target.id] = set(
                        pathway for pathway in universal_model.metabolite_pathway_map[target.id])
        pathways_to_keep = set(pathways_to_keep) - pathways_to_ignore
        metabolite_pathways_map = {metabolite: pathways - pathways_to_ignore for metabolite, pathways in
                                   metabolite_pathways_map.items() if pathways not in pathways_to_ignore}
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
        print('Number of reactions in temporary universal model:',
              len(self.universal_model_compartmentalized.reactions))
        write_sbml_model(self.universal_model_compartmentalized,
                         os.path.join(folder_path, 'temporary_universal_model.xml'))

    def identify_additional_seeds(self, combinet, targets):
        """
        Identify additional seeds to add to the model.
        Parameters
        ----------
        combinet : TermSet
            The combinet model, Combined model with custom universal model and draft model.
        targets : cobra.Model
            The targets model.


        Returns
        -------
        """
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
                                self.additional_seeds.add(
                                    self.universal_model.metabolites.get_by_id(potential_cofactor))
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
        """
        Add demands to the model.
        Returns
        -------

        """
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
            self.original_model.add_reactions(
                self.universal_model.reactions.get_by_id(rxn.id) for rxn in list(self.positive_solutions)[0])

    def find_deadends(self, model=None):
        """
        Find dead-end metabolites in the model.
        Parameters
        ----------
        model : cobra.Model, optional
            The model to use. If None, the original model will be used.

        Returns
        -------
        dead_ends : list
            List of dead-end metabolites.
        """

        if model is None:
            model = self.original_model

        def is_only_product(metabolite: cobra.Metabolite, reaction: cobra.Reaction) -> bool:
            """Determine if a metabolite is only a product of a reaction.

            Parameters
            ----------
            metabolite : cobra.Metabolite
                The metabolite to check.
            reaction : cobra.Reaction
                The reaction to check.

            Returns
            -------
            bool
                True if the metabolite is only a product of the reaction, False otherwise.
            """

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
        self.required_additional_seeds_and_demands = self.additional_seeds.union(
            set([met.id for demand in self.original_model.demands for met in demand.metabolites]))
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
            self.required_additional_seeds_and_demands = self.required_additional_seeds_and_demands - set(
                [met.id for met in sink.metabolites])

    def add_reactions_to_model_v2(self):
        """
        Add reactions to a model based on the obtained completion.
        For each optimal completion, it adds the reactions. If a positive solution is obtained (fba sol) it adds the solution to a list.
        Otherwise, it adds demands to the model and tries to optimize again. If a positive solution is obtained, it adds the solution to a list.
        In the end, if at least one positive solution is obtained, the first one is selected as the final model. Otherwise, the minimal completion
        is added to the model (even without producing biomass).
        Returns
        -------

        """
        positive_solutions = set()
        already_in_original_model = {reaction.id for reaction in self.original_model.reactions}
        seeds = read_sbml_model(self.seeds_path)
        mets_in_medium = self.get_mets_in_medium(self.original_model)

        for m in self.optimal_completions:
            completion = self.get_completion(m)
            self.all_completions.append(list(completion))
            temp_model = self.original_model.copy()
            to_add, already_in_original_model = self.add_reactions_to_temp_model(completion, already_in_original_model,
                                                                                 temp_model)

            self.create_sinks_for_metabolites(temp_model, seeds, mets_in_medium)

            self.check_for_positive_solution(temp_model, to_add, positive_solutions)

        self.handle_final_model_selection(positive_solutions, seeds, mets_in_medium, already_in_original_model)

    def get_mets_in_medium(self, model):
        """
        Get metabolites in the medium.
        Parameters
        ----------
        model : cobra.Model
            The model to use.

        Returns
        -------

        """
        return {exchange.reactants[0].id for exchange in model.exchanges if exchange.lower_bound < 0}

    def get_completion(self, m):
        """
        Get the completion from the gap-filling results.
        Parameters
        ----------
        m : dict
            The gap-filling results.

        Returns
        -------

        """
        return set(a[0] for pred in m if pred == 'xreaction' for a in m[pred])

    def add_reactions_to_temp_model(self, completion, already_in_original_model, temp_model):
        """
        Add reactions to a temporary model based on a given set of reactions (completion).
        It also updates the set of reactions already present in the original model.

        Parameters:
        - completion (set): A set of reactions to be added to the model.
        - already_in_original_model (set): A set of reaction IDs already in the original model.
        - temp_model (cobra.Model): The temporary model to which reactions will be added.

        Returns:
        - tuple: A tuple containing the list of added reactions and the updated set of reactions already in the original model.
        """
        to_add = []
        for reaction in completion:
            save_reaction = self.universal_model.reactions.get_by_id(reaction.replace("R_", "")).copy()
            if save_reaction.id in already_in_original_model:
                save_reaction.id = save_reaction.id + "_universal"
            to_add.append(save_reaction)

        temp_model.add_reactions(to_add)
        # Ensure already_in_original_model is treated as a set before updating
        already_in_original_model = set(already_in_original_model)
        already_in_original_model.update({r.id for r in to_add})
        return to_add, already_in_original_model

    def create_sinks_for_metabolites(self, model, seeds, mets_in_medium):
        """
        Create sinks for metabolites in the medium.
        Parameters
        ----------
        model : cobra.Model
            The model to use.
        seeds : cobra.Model
            The model containing seed metabolites.
        mets_in_medium : set
            A set of metabolite IDs that are present in the medium.

        Returns
        -------

        """

        model_ids = [e.id for e in model.metabolites]
        for seed in seeds.metabolites:
            if seed.id in model_ids and not {seed.id}.issubset(mets_in_medium):
                model.create_sink(seed.id)  # type: ignore

    def check_for_positive_solution(self, temp_model, to_add, positive_solutions):
        """
        Check if the temporary model has a positive solution. If so, add it to the set of positive solutions.
        Parameters
        ----------
        temp_model : cobra.Model
            The temporary model to check.
        to_add : list
            List of reactions that were added to the model.
        positive_solutions : set
            A set to store models with positive solutions.

        Returns
        -------

        """
        sol = temp_model.slim_optimize()
        if round(sol, 5) > 0:
            positive_solutions = set(positive_solutions)
            positive_solutions.add((temp_model, tuple([reaction.id for reaction in to_add])))
        else:
            self.add_demands_and_reoptimize(temp_model, to_add, positive_solutions)

    def add_demands_and_reoptimize(self, temp_model, to_add, positive_solutions):
        """
        Add demand reactions to the temporary model and reoptimize. If a positive solution is found,
        add it to the set of positive solutions.
        Parameters
        ----------
        temp_model : cobra.Model
            The temporary model to which demands are added and reoptimized.
        to_add : list
            List of reactions that were added to the model.
        positive_solutions : set
            A set to store models with positive solutions.

        Returns
        -------

        """
        temp_model, sol, demands = self.add_demands_v2(temp_model)  # type: ignore
        if round(sol, 5) > 0:
            to_add += demands
            # Ensure positive_solutions is treated as a set before using the add method
            positive_solutions = set(positive_solutions)
            positive_solutions.add((temp_model, tuple([reaction.id for reaction in to_add])))

    def handle_final_model_selection(self, positive_solutions, seeds, mets_in_medium, already_in_original_model):
        """
        Select the final model based on the available positive solutions. If no positive solution exists,
        use the minimal completion.
        Parameters
        ----------
        positive_solutions : set
            A set containing models with positive solutions.
        seeds : cobra.Model
            The model containing seed metabolites.
        mets_in_medium : set
            A set of metabolite IDs that are present in the medium.
        already_in_original_model : set
            A set of reaction IDs already in the original model.

        Returns
        -------

        """

        if positive_solutions:
            self.original_model = list(positive_solutions)[0][0]
            print(self.original_model.slim_optimize())
        else:
            self.handle_no_positive_solution(seeds, mets_in_medium, already_in_original_model)

    def handle_no_positive_solution(self, seeds, mets_in_medium, already_in_original_model):
        """

        Handle the case when no positive solution is found by adding minimal completion to the model.

        Parameters
        ----------
        seeds : cobra.Model
            The model containing seed metabolites.
        mets_in_medium : set
            A set of metabolite IDs that are present in the medium.
        already_in_original_model : set
            A set of reaction IDs already in the original model.

        Returns
        -------

        """
        to_add = []
        for reaction in self.minimal_completion:
            save_reaction = self.universal_model.reactions.get_by_id(reaction.replace("R_", "")).copy()
            if {save_reaction.id}.issubset(already_in_original_model):
                save_reaction.id = save_reaction.id + "_universal"
            to_add.append(save_reaction)
        self.create_sinks_for_metabolites(self.original_model, seeds, mets_in_medium)
        self.original_model.add_reactions(to_add)
