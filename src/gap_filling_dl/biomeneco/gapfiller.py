import os
import time
import cobra
from meneco import run_meneco


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

        return gap_filler

    # def evaluate_results(self, verbose=False):
    #     """
    #     Evaluate the results of the gap-filling algorithm.
    #
    #     Parameters
    #     ----------
    #     verbose : bool, optional
    #         Whether to print out detailed information about the results.
    #
    #     Returns
    #     -------
    #     dict
    #         A dictionary containing the evaluation metrics.
    #     """
    #     results = {}
    #
    #     # Get the number of seed reactions
    #     seeds = cobra.io.read_sbml_model(self.seeds_path)
    #     results["Seed reactions"] = len(seeds.reactions)
    #
    #     # Get the number of target reactions
    #     targets = cobra.io.read_sbml_model(self.targets_path)
    #     results["Target reactions"] = len(targets.reactions)
    #
    #     # Get the number of reactions in the draft network
    #     draft = cobra.io.read_sbml_model(self.results_meneco["Draft network file"])
    #     results["Draft reactions"] = len(draft.reactions)
    #
    #     # Get the number of gap-filled reactions
    #     gf = cobra.io.read_sbml_model(self.results_meneco["Draft network file"])
    #     gf.remove_reactions([r.id for r in draft.reactions])
    #     results["Gap-filled reactions"] = len(gf.reactions)
    #
    #     # Get the number of unproducible targets
    #     results["Unproducible targets"] = len(self.results_meneco["Unproducible targets"])
    #
    #     if verbose:
    #         print("Evaluation results:")
    #         for key, value in results.items():
    #             print(f"\t{key}: {value}")
    #
    #     return results

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

