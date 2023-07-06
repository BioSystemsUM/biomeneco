import os
import time
import cobra
from meneco.query import get_minimal_completion_size
from meneco.meneco import run_meneco
from meneco.meneco import sbml
import copy

from gap_filling_dl.utils import build_temporary_universal_model


class GapFiller:
    def __init__(self, model_path, universal_model_path):

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

            # Call the get_minimal_completion_size function this order: draftnet, repairnet, seeds, targets)
            self.minimal_completion = get_minimal_completion_size(draftnet, repairnet, seeds, targets)
            time_end = time.time()
            self.gftime = time_end - time_start

            self.minimal_set = set(
                atom[0] for pred in self.minimal_completion if pred == 'xreaction'
                for atom in self.minimal_completion[pred])

            print("Time taken: ", self.gftime)

        return self.report(write_to_file=write_to_file, removed_reactions=removed_reactions,
                           objective_function_id=objective_function_id)

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
                    **kwargs):
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
        t = os.listdir(folder_path)
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

        if temporary_universal_model:
            if not objective_function_id:
                raise ValueError(
                    "Objective function ID must be specified for the creation of a temporary universal model.")

            build_temporary_universal_model(gap_filler, folder_path, related_pathways=True)

            if os.path.isfile(os.path.join(folder_path, 'temporary_universal_model.xml')):
                print('Temporary universal model file successfully created.')
                gap_filler.universal_model_path = os.path.join(folder_path, 'temporary_universal_model.xml')
                print('Temporary universal model file path:', gap_filler.universal_model_path)
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

    def report(self, write_to_file=False, removed_reactions: list = None, objective_function_id: str = None, ):

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

                #universal_model.reactions.get_by_id('

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
