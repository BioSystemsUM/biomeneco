def fix_terms(input_file, output_file):
    with open(input_file, 'r') as f:
        data = f.read()

    # Replace 'extr' with 'out'
    data = data.replace('extr', 'out')

    # Replace 'cytop' with 'in'
    data = data.replace('cytop', 'in')

    with open(output_file, 'w') as f:
        f.write(data)


# path = '/Users/josediogomoura/gap_filling_dl/tests/data/6gaps/6gaps_model.xml'
# fix_terms(path, path)
#
# path_seeds = '/Users/josediogomoura/gap_filling_dl/tests/data/6gaps/6gaps_seeds.xml'
# fix_terms(path_seeds, path_seeds)
#
# path_targets = '/Users/josediogomoura/gap_filling_dl/tests/data/6gaps/6gaps_targets.xml'
# fix_terms(path_targets, path_targets)


