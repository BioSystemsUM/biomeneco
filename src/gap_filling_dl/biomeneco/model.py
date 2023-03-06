from typing import Tuple, List
import cobra
from bioiso import BioISO
from bioiso import load, set_solver
from src.gap_filling_dl.write_sbml.metabolites_to_sbml import write_metabolites_to_sbml


class Model(cobra.Model):
    def __init__(self, model_path, *args, **kwargs):
        self.model_path = model_path
        self.initialize_model()
        super().__init__(*args, **kwargs)

    def initialize_model(self):
        """
        Initialize the model.

        """
        self.model = load(self.model_path)

    def __str__(self):
        """
        Print the model summary.
        Returns a string with the model summary.
        -------

        """

        return self.model.summary()

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

    def identify_targets(self, objective_function: str = None,
                         objective='maximize', solver='cplex') -> List[Tuple[str, str]]:
        """
    Identify the targets (external metabolites that leave the model) of a model.

    Parameters
    ----------
    objective_function: str
        The name of the objective function to use. Default is 'Biomass_C3_cytop'.
    objective: str
        The type of optimization objective to use. Default is 'maximize'.
    solver: str
        The solver to use for optimization. Default is 'cplex'.

    Returns
    -------
    List of tuples representing identified targets, where each tuple is of the form (metabolite_id, compartment).
    """

        if objective_function is None:  # 'Biomass_C3_cytop'
            raise ValueError("Objective function cannot be None.")

        if objective is None:
            raise ValueError("Objective type cannot be None.")

        if solver is None:
            raise ValueError("Solver cannot be None.")

        model = self.model
        # print(model.summary())
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

    def to_sbml(self, file_name, save_path, seeds=True, targets=True):
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
            write_metabolites_to_sbml(file_name, save_path, self.seeds)
        if targets:
            write_metabolites_to_sbml(file_name, save_path, self.targets)







