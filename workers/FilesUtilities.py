from zipfile import ZipFile
import os
import subprocess


def read_conf_file(path):
    '''
    reading the configuration file
    It must be assembled in the following manner: key = configuration value
    It will return a dictionary in the following format: {workers : 3, limit : 100}
    '''
    res = {}

    with open(path) as file:
        lines = file.readlines()
        for line in lines:
            if '=' in line:
                prop = line.replace(" ", "").replace("\n", "")
                l_str = prop.split('=')
                res[l_str[0]] = l_str[1]

    return res


def read_parameters_file(processingPath, resultsPath, SubmissionParameters_path):
    value_parameters = ["db", "outfmt", "strand", "unal", "max-target-seqs", "top", "max-hsps", "range-culling", \
                        "evalue", "min-score", "id", "query-cover", "subject-cover", "block-size", "index-chunks", "gapopen", \
                        "gapextend", "frameshift", "long-reads", "matrix", "comp-based-stats", "masking", "no-self-hits", \
                        "taxonlist", "taxon-exclude", "seqidlist", "algo", "bin", "min-orf", "freq-sd", "id2", "xdrop", "shapes", \
                        "shape-mask", "ext", "culling-overlap", "taxon-k", "range-cover"]
    in_file_parameters = ["custom-matrix"]
    non_value_parameters = ["mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive", "skip-missing-seqids", \
                            "xml-blord-format"]

    taxonomy_parameters = ["db", "outfmt", "max-target-seqs", "top", "max-hsps", "gapopen", "gapextend", "matrix", "shapes", "shape-mask", \
                           "custom-matrix", "mid-sensitive", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"]

    parameters = {}

    with open(SubmissionParameters_path) as file:
        lines = file.readlines()
        for line in lines:
            if '=' in line:
                prop = line.replace("\n", "")
                l_str = prop.split(" = ")
                parameter = l_str[0].replace(" ", "")
                value = l_str[1]

                if parameter in non_value_parameters:
                    parameters[parameter] = ""

                elif parameter in in_file_parameters:
                    file_path = processingPath + value
                    if os.path.exists(file_path):
                        parameters[parameter] = value
                    elif os.path.exists(resultsPath + value):
                        parameters[parameter] = value
                    else:
                        parameters[parameter] = None

                elif parameter in value_parameters:
                    parameters[parameter] = value

                elif parameter == "taxonomyRankID":
                    filename = processingPath + value

                    if os.path.exists(filename):
                        ids = read_taxids(filename)
                        parameters["taxonlist"] = ids
                    else:
                        subprocess.call("get_species_taxids.sh -t " + value + " > " + filename, shell=True)
                        if os.path.exists(filename):
                            ids = read_taxids(filename)
                            parameters["taxonlist"] = ids

                else:
                    parameters[parameter] = None

    parameters_new = {}

    if "db" in parameters.keys():
        if parameters["db"] == "UniRef50":
            for parameter in parameters.keys():
                if parameter in taxonomy_parameters:
                    parameters_new[parameter] = parameters[parameter]

            parameters = parameters_new

    parameters["threads"] = "3"

    return parameters


def read_taxids(filepath):
    ids = ""
    with open(filepath, "r") as file:
        for line in file.readlines():
            ids = ids + "," + line.replace("\n", "")
    return ids


def read_workers_conf_file(path):
    res = []

    with open(path) as file:
        lines = file.readlines()
        for line in lines:
            if "#" not in line[0]:
                str = line.replace(" ", "").replace("\n", "")
                res.append(str)

    return res


def compressFiles(path, savePath):
    file_paths = get_all_file_paths(path)

    zipLocation = savePath

    with ZipFile(zipLocation, 'w') as zip:
        # writing each file one by one
        for file in file_paths:
            new_file = file.split("/")[-1]
            zip.write(file, arcname=new_file)

    return zipLocation


def get_all_file_paths(directory):
    # initializing empty file paths list
    file_paths = []

    # crawling through directory and subdirectories
    for root, directories, files in os.walk(directory):

        for filename in files:
            # join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)

    # returning all file paths
    return file_paths