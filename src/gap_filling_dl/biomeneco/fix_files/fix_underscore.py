import tempfile
import sys


def fix_file(input_file):
    with open(input_file, "r") as f:
        data = f.read()

    # Replace "__" with "_"
    data = data.replace("__", "_")

    with open(input_file, "w") as f:
        f.write(data)


# try with a file in a path:
#fix_file("/Users/josediogomoura/gap_filling_dl/tests/performance_tests/data/kegg_universal_model.xml")
fix_file("/Users/josediogomoura/gap_filling_dl/tests/data/universal_model_kegg.xml")

