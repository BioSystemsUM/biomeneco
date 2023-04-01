import unittest
import os
import tempfile
from src.gap_filling_dl.biomeneco.fix_files.fix_terms import fix_terms


class TestFixTerms(unittest.TestCase):

    def setUp(self):
        self.input_file = tempfile.NamedTemporaryFile(delete=False).name
        self.output_file = tempfile.NamedTemporaryFile(delete=False).name

    def test_fix_terms(self):
        # Create test input file
        with open(self.input_file, 'w') as f:
            f.write('extracellular cytoplasmic')

        # Call function to fix terms
        fix_terms(self.input_file, self.output_file)

        # Check if output file exists
        self.assertTrue(os.path.exists(self.output_file))

        # Check if output file contains correct replacements
        with open(self.output_file, 'r') as f:
            output_data = f.read()
        self.assertEqual(output_data, 'outacellular inlasmic')

    def tearDown(self):
        # Remove test files
        os.remove(self.input_file)
        os.remove(self.output_file)


if __name__ == '__main__':
    unittest.main()
