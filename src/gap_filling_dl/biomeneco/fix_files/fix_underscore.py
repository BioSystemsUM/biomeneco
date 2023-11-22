import tempfile
import sys


def fix_file(input_file):
    with open(input_file, "r") as f:
        data = f.read()

    # Replace "__" with "_"
    data = data.replace("__", "_")

    with open(input_file, "w") as f:
        f.write(data)

