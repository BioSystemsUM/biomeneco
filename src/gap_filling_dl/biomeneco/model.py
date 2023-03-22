import os
from typing import Tuple, List
import cobra
from bioiso import BioISO
from bioiso import load, set_solver
from src.gap_filling_dl.write_sbml.metabolites_to_sbml import write_metabolites_to_sbml


class Model(cobra.Model):
    def __init__(self, model: cobra.Model, objective_function_id, *args, **kwargs):
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
                if not specific_biomass_components and biomass_components[biomass_component].get("identifier") not in targets:
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

        if seeds:  # if the user wants to save the seeds.xml
            seeds_file_name = os.path.splitext(file_name)[0] + "_seeds.xml"
            write_metabolites_to_sbml(seeds_file_name, save_path, self.seeds)
        if targets:  # if the user wants to save the targets
            targets_file_name = os.path.splitext(file_name)[0] + "_targets.xml"
            write_metabolites_to_sbml(targets_file_name, save_path, self.targets)
