import unittest
import os
from tempfile import TemporaryDirectory
from gap_filling_dl.biomeneco.fix_files.fix_underscore import fix_file


class TestFixFile(unittest.TestCase):
    def test_fix_file(self):
        # Create a temporary directory and file to test on
        with TemporaryDirectory() as tmpdir:
            input_file = os.path.join(tmpdir, "test.txt")
            with open(input_file, "w") as f:
                f.write("hello__world")

            # Call the function being tested
            fix_file(input_file)

            # Read the contents of the modified file
            with open(input_file, "r") as f:
                result = f.read()

            # Check that the function made the expected changes
            expected = "hello_world"
            self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
