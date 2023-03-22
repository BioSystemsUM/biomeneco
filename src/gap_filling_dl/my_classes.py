import os
from typing import Tuple, List
import time
import cobra
from cobra.io import write_sbml_model
from cobra import Metabolite
from bioiso import BioISO
from bioiso import load, set_solver
from meneco import run_meneco


def write_metabolites_to_sbml(file_name, save_path, metabolites):
    """
    Write a list of metabolites to an SBML file.

    Parameters
    ----------
    file_name: str
        The name of the SBML file to write to.
    save_path: str
        The path to save the SBML file to.
    metabolites: List[Tuple[str, str]]
        A list of metabolites to write to the SBML file.

    """

    model = cobra.Model()
    metabolite_cobra_metabolite_objects = []
    for metabolite_id, compartment in metabolites:
        metabolite = Metabolite(id=metabolite_id, compartment=compartment)
        metabolite_cobra_metabolite_objects.append(metabolite)

    model.add_metabolites(metabolite_cobra_metabolite_objects)
    write_sbml_model(model, os.path.join(save_path, file_name))
    

class Model(cobra.Model):
    def __init__(self, model_path, *args, **kwargs):
        self.model_path = model_path
        self.initialize_model()
        super().__init__(*args, **kwargs)

    def initialize_model(self):
        self.model = load(self.model_path)

    def __str__(self):
        return self.model.summary()

    @property
    def seeds(self) -> List[Tuple[str, str]]:
        return self.identify_seeds()

    @property
    def targets(self) -> List[Tuple[str, str]]:
        return self.identify_targets()

    # method to identify the seeds of a model
    def identify_seeds(self) -> List[Tuple[str, str]]:
        """
        Identify the seeds (external metabolites that enter the model) of a model.

        Parameters
        ----------
        model_path: str
            The path to the model to identify the seeds of.
        """

        model = self.model
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
                    if reaction.upper_bound > 0 and product.id not in total_seeds_ids:
                        seeds_ids.append(product.id)
                        seeds.append((product.id, product.compartment))

            if found_out and found_inside:
                total_seeds += seeds
                total_seeds_ids += seeds_ids

        return total_seeds

    def identify_targets(self, objective_function="Biomass_C3_cytop",
                         objective="maximize", solver='cplex') -> List[Tuple[str, str]]:
        """
        Identify the targets (external metabolites that leave the model) of a model.

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
        """
        model = self.model
        #print(model.summary())
        set_solver(model, solver)
        model.objective = objective_function

        bio = BioISO(objective_function, model, objective)
        bio.run(2, False)

        results = bio.get_tree()
        biomass_components = results["M_root_M_root_M_root_product"]["next"]

        targets = []

        for biomass_component in biomass_components:

            biomass_component_role = biomass_components[biomass_component].get("role")
            biomass_component_analysis = biomass_components[biomass_component].get("analysis")
            if biomass_component_role == "Reactant" and not biomass_component_analysis:
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

    def to_sbml(self, file_name, save_path, seeds=None, targets=None):
        """
        Write a model to an SBML file.

        Parameters
        ----------
        file_name: str
            The name of the SBML file to write to.
        save_path: str
            The path to save the SBML file to.
        seeds: list[Tuple[str, str]]
            A list of tuples of seed metabolite IDs and compartments.
        targets:
            A list of tuples of target metabolite IDs and compartments.
        """

        if seeds:
            write_metabolites_to_sbml(self.model_path, file_name, save_path, seeds)
        if targets:
            write_metabolites_to_sbml(self.model_path, file_name, save_path, targets)


class GapFiller:
    def __init__(self, model_path, universal_model_path):

        self.results_meneco = None
        self.universal_model_path = universal_model_path
        self.model_path = model_path
        self._seeds_path = None
        self._targets_path = None

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
    def universal_model(self):
        return cobra.io.read_sbml_model(self.universal_model_path)

    def run(self, **kwargs):
        """
        Run the gap-filling algorithm.

        Parameters
        ----------
        kwargs:
            Keyword arguments to pass to the gap-filling algorithm.

        """
        time_start = time.time()
        self.run_meneco(**kwargs)

        time_end = time.time()
        print("Time taken: ", time_end - time_start)

        return self.results_meneco

    def run_meneco(self, enumeration=True, json_output=True):
        """
        Run the gap-filling algorithm.

        Parameters
        ----------
        kwargs:
            Keyword arguments to pass to the gap-filling algorithm.

        Returns
        -------
        results: dict
        """
        self.results_meneco = run_meneco(self.model_path, self.universal_model_path, self.seeds_path, self.targets_path,
                                         enumeration=enumeration, json_output=json_output)
        return self.results_meneco

    # @classmethod
    # def from_folder(cls, folder_path):
    #     """
    #     Create a GapFiller object from a folder.
    #
    #     Parameters
    #     ----------
    #     folder_path: str
    #         The path to the folder to create a GapFiller object from.
    #     """
    #     gap_filler = cls(folder_path + "/model.xml", folder_path + "/universal_model.xml")
    #     gap_filler.seeds_path = folder_path + "/seeds.xml"
    #     gap_filler.targets_path = folder_path + "/targets.xml"
    #     return gap_filler

    @classmethod
    def from_folder(cls, folder_path):
        """
        Create a GapFiller object from a folder.

        Parameters
        ----------
        folder_path: str
            The path to the folder to create a GapFiller object from.
        """

        model_file = None
        seeds_file = None
        targets_file = None
        universal_model_file = None

        for file in os.listdir(folder_path):
            if file.endswith(".xml"):
                if "model" in file:
                    model_file = file
                elif "seeds" in file:
                    seeds_file = file
                elif "targets" in file:
                    targets_file = file
                elif "universal_model" in file:
                    universal_model_file = file

        if not model_file:
            raise FileNotFoundError("No model file found in folder.")
        if not seeds_file:
            raise FileNotFoundError("No seeds file found in folder.")
        if not targets_file:
            raise FileNotFoundError("No targets file found in folder.")
        if not universal_model_file:
            raise FileNotFoundError("No universal model file found in folder.")

        gap_filler = cls(folder_path + "/" + model_file, folder_path + "/" + universal_model_file)
        gap_filler.seeds_path = folder_path + "/" + seeds_file
        gap_filler.targets_path = folder_path + "/" + targets_file
        return gap_filler

    def evaluate_results(self, verbose=False):
        """
        Evaluate the results of the gap-filling algorithm.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print out detailed information about the results.

        Returns
        -------
        dict
            A dictionary containing the evaluation metrics.
        """
        results = {}

        # Get the number of seed reactions
        seeds = cobra.io.read_sbml_model(self.seeds_path)
        results["Seed reactions"] = len(seeds.reactions)

        # Get the number of target reactions
        targets = cobra.io.read_sbml_model(self.targets_path)
        results["Target reactions"] = len(targets.reactions)

        # Get the number of reactions in the draft network
        draft = cobra.io.read_sbml_model(self.results_meneco["Draft network file"])
        results["Draft reactions"] = len(draft.reactions)

        # Get the number of gap-filled reactions
        gf = cobra.io.read_sbml_model(self.results_meneco["Draft network file"])
        gf.remove_reactions([r.id for r in draft.reactions])
        results["Gap-filled reactions"] = len(gf.reactions)

        # Get the number of unproducible targets
        results["Unproducible targets"] = len(self.results_meneco["Unproducible targets"])

        if verbose:
            print("Evaluation results:")
            for key, value in results.items():
                print(f"\t{key}: {value}")

        return results





# fazer from folder: Ok
# fazer tests ok~~~
# alterar id das reações e metabolitos do universal model ok
# docstrings, tipagem: quase tudo
# dividir classes em ficheiros diferentes: Ok















