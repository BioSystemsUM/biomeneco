import unittest
import os
import tempfile
import cobra
from src.gap_filling_dl.write_sbml.metabolites_to_sbml import write_metabolites_to_sbml


class TestWriteMetabolitesToSBML(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.save_path = self.tempdir
        self.file_name = "test_model.xml"
        self.model = cobra.Model()
        self.model.add_metabolites(
            [
                cobra.Metabolite("M_C00054_cytop", compartment="cytop"),
                cobra.Metabolite("M_C00008_cytop", compartment="cytop"),
                cobra.Metabolite("M_C00750_extr_b", compartment="extr_b"),
            ]
        )
        self.model.add_reactions(
            [
                cobra.Reaction("R1"),
                cobra.Reaction("R2"),
                cobra.Reaction("R3"),
            ]
        )
        self.model.reactions.get_by_id("R1").add_metabolites(
            {
                self.model.metabolites.get_by_id("M_C00054_cytop"): -1.0,
                self.model.metabolites.get_by_id("M_C00008_cytop"): 1.0,
            }
        )
        self.model.reactions.get_by_id("R2").add_metabolites(
            {
                self.model.metabolites.get_by_id("M_C00008_cytop"): -1.0,
                self.model.metabolites.get_by_id("M_C00750_extr_b"): 1.0,
            }
        )
        self.model.reactions.get_by_id("R3").add_metabolites(
            {
                self.model.metabolites.get_by_id("M_C00054_cytop"): -1.0,
                self.model.metabolites.get_by_id("M_C00750_extr_b"): 1.0,
            }
        )
        self.metabolites = [
            ("M_C00054_cytop", "cytop"),
            ("M_C00008_cytop", "cytop"),
            ("M_C00750_extr_b", "extr_b"),
        ]

    def tearDown(self):
        for file in os.listdir(self.tempdir):
            os.remove(os.path.join(self.tempdir, file))
        os.rmdir(self.tempdir)

    def test_write_metabolites_to_sbml(self):
        write_metabolites_to_sbml(self.file_name, self.save_path, self.metabolites)
        self.assertTrue(os.path.exists(os.path.join(self.save_path, self.file_name)))

        # Check if the model is not empty
        self.assertTrue(len(self.model.reactions) > 0)

        # Check if the expected metabolites are in the model
        metabolite_ids = set([met.id for met in self.model.metabolites])
        expected_ids = set([m[0] for m in self.metabolites])
        self.assertEqual(metabolite_ids, expected_ids)


if __name__ == "__main__":
    unittest.main()












# def test_write_metabolites_to_sbml():
#     # Create temporary directory for testing
#     tempdir = tempfile.mkdtemp()
#     # Set the save path to the temporary directory
#     save_path = tempdir
#     # Set the file name
#     file_name = "test_model.xml"
#     # Create a list of metabolites
#     metabolites = [("glc_D_c", "c"), ("atp_c", "c")]
#     # Call the function
#     write_metabolites_to_sbml(file_name, save_path, metabolites)
#     # Check that the file was created
#     assert os.path.exists(os.path.join(save_path, file_name))
#     # Load the created file with COBRApy
#     model = cobra.io.read_sbml_model(os.path.join(save_path, file_name))
#     # Check that the model has the expected metabolites
#     metabolite_ids = set([met.id for met in model.metabolites])
#     expected_ids = set([m[0] for m in metabolites])
#     assert metabolite_ids == expected_ids
#     # Clean up the temporary directory
#     os.rmdir(tempdir)
#
#
# if __name__ == '__main__':
#     unittest.main()
