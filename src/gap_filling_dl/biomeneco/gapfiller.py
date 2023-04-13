import os
import time
import cobra
from meneco import run_meneco
from cobra.io import read_sbml_model, write_sbml_model
from src.gap_filling_dl.biomeneco.model import Model


class GapFiller:
    def __init__(self, model_path, universal_model_path):

        self.cobra_filled_model = None
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
    def from_folder(cls, folder_path, temporary_universal_model=False, objective_function_id=None):
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
        # print(os.listdir(folder_path))
        # print(os.listdir(folder_path))
        for file in os.listdir(folder_path):
            if file.endswith(".xml"):
                if "model" in file and "universal" not in file:
                    model_file = file
                elif "seeds" in file:
                    seeds_file = file
                elif "targets" in file:
                    targets_file = file
                elif "universal_model" in file:
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

        # create the GapFiller object
        gap_filler = cls(model_path=os.path.join(folder_path, model_file),
                         universal_model_path=os.path.join(folder_path, universal_model_file))
        gap_filler.seeds_path = os.path.join(folder_path, seeds_file)
        gap_filler.targets_path = os.path.join(folder_path, targets_file)

        # create a temporary universal model if necessary

        if temporary_universal_model:
            pathways_to_ignore = []
            pathways_to_keep = []

            # create a Model object
            my_model = Model(model=gap_filler.model_path, objective_function_id=objective_function_id)

            # get the pathways to ignore
            for target in gap_filler.targets_path:
                pathways_to_keep += [pathway for pathway in my_model.metabolite_pathway_map[target[0]]]
            pathways_to_keep = list(set(pathways_to_keep))
            print('Pathways to keep are:', pathways_to_keep)
            universal_model = cobra.io.read_sbml_model(gap_filler.universal_model_path)
            print('Number of reactions in universal model:', len(universal_model.reactions))
            to_keep = set()

            # remove the reactions that are not in the pathways to keep
            for pathway in universal_model.groups:
                if pathway.name in pathways_to_keep:
                    to_keep.update(reaction for reaction in pathway.members)
            to_remove = set(universal_model.reactions) - to_keep
            universal_model.remove_reactions(list(to_remove), remove_orphans=True)
            print('Number of reactions in temporary universal model:', len(universal_model.reactions))

            # overwrite the universal model file? or create a new one?
            write_sbml_model(universal_model, os.path.join(folder_path, 'universal_model.xml'))
            gap_filler.universal_model_path = os.path.join(folder_path, 'universal_model.xml')

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

        # Get the number of unproducible reactions
        results["Unproducible reactions"] = len(self.results_meneco["Unproducible targets"])

        if verbose:
            print("Evaluation results:")
            for key, value in results.items():
                print(f"\t{key}: {value}")

        return results

    def add_reactions_to_model(self, reactions: list[str]):
        """
        Add reactions to a model and run the cobra gap-filling algorithm.

        Parameters
        ----------
        reactions : list
            A list of reactions to add to the model.
        """
        model = cobra.io.read_sbml_model(self.model_path)
        uni_model = cobra.io.read_sbml_model(self.universal_model_path)

        # filter the input reactions to remove reactions not found in the universal model
        reactions_to_add = [reaction for reaction in reactions if reaction.replace("R_", "") in uni_model.reactions]

        # add the filtered reactions to the model
        model.add_reactions(
            [uni_model.reactions.get_by_id(reaction.replace("R_", "")) for reaction in reactions_to_add])

        # run the cobra gap-filling algorithm
        self.cobra_filled_model = cobra.flux_analysis.gapfilling.gapfill(model, uni_model, demand_reactions=True,
                                                                         exchange_reactions=True)

        return self.cobra_filled_model

    def compare_reactions(self):
        """
        Compare the reactions in the initial model and the filled model.
        Parameters
        ----------

        Returns
        -------

        """

        # Get the reactions in the initial model
        initial_reactions = cobra.io.read_sbml_model(self.model_path).reactions

        # Get the reactions in the filled model
        filled_model_reactions = self.cobra_filled_model.reactions

        # Get the IDs of the reactions in each model
        initial_reaction_ids = set(r.id for r in initial_reactions)
        filled_model_reaction_ids = set(r.id for r in filled_model_reactions)

        # Find the reactions that were added to the model during the gap-filling process
        added_reactions = list(filled_model_reaction_ids - initial_reaction_ids)

        return added_reactions
