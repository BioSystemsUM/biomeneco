def fix_terms(input_file, output_file):
    with open(input_file, 'r') as f:
        data = f.read()

    # Replace 'extr' with 'out'
    data = data.replace('extr', 'out')

    # Replace 'cytop' with 'in'
    data = data.replace('cytop', 'in')

    with open(output_file, 'w') as f:
        f.write(data)
