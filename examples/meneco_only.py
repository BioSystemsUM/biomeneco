import logging
import sys
import time
from clyngor.as_pyasp import TermSet, Atom
from cobra.io import read_sbml_model, write_sbml_model
from gap_filling_dl import write_metabolites_to_sbml

from gap_filling_dl.biomeneco.model import Model
from meneco.meneco import query, sbml, extract_xreactions
import warnings

logging.getLogger("cobra").setLevel(logging.WARNING)
warnings.filterwarnings("ignore")
logging.getLogger("cobra").setLevel(logging.CRITICAL)
logging.getLogger("meneco").setLevel(logging.CRITICAL)
import platform


if platform.system() == 'Windows':
    prefix = "../"
else:
    prefix = "/home/"


def run_meneco(paths_map):
    time_start = time.time()
    draftnet = sbml.readSBMLnetwork(paths_map['model'], "draft")

    seeds = sbml.readSBMLseeds(paths_map['seeds'])

    targets = sbml.readSBMLtargets(paths_map['targets'])

    model = query.get_unproducible(draftnet, targets, seeds)
    unproducible = set([a[0] for pred in model if pred == 'unproducible_target' for a in model[pred]])
    print('{0} unproducible targets:\n\t{1}\n'.format(len(unproducible), '\n\t'.join(unproducible)))

    repairnet = sbml.readSBMLnetwork(paths_map['repairnet'], 'repairnet')

    combinet = draftnet
    combinet = TermSet(combinet.union(repairnet))

    additional_seeds = sbml.readSBMLseeds(paths_map['additional_seeds'])
    model = query.get_unproducible(combinet, targets, additional_seeds)
    never_producible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
    print('{0} unreconstructable targets:\n\t{1}\n'.format(len(never_producible), '\n\t'.join(never_producible)))

    reconstructable_targets = set(list(unproducible.difference(never_producible)))
    print('{0} reconstructable targets:\n\t{1}\n'.format(len(reconstructable_targets), '\n\t'.join(reconstructable_targets)))
    #

    targets = TermSet(Atom('target("' + t + '")') for t in reconstructable_targets)
    optimum = extract_xreactions(query.get_minimal_completion_size(draftnet, repairnet, additional_seeds, targets), False)

    print('minimal size =', optimum)

    print('One minimal completion of size {0}:\n\t{1}\n'.format(
        optimum, '\n\t'.join(optimum)))

    optimal_completions = query.get_optimal_completions(draftnet, repairnet, seeds, targets, len(optimum))
    exec_time = time.time() - time_start

    print('execution time = {0}'.format(exec_time))

    cobra_model = read_sbml_model(paths_map['model'])
    universal_model = read_sbml_model(paths_map['repairnet'])

    cobra_model.add_reactions([universal_model.reactions.get_by_id(reaction.replace("R_", "")) for reaction in optimum])
    sol = cobra_model.slim_optimize()
    # all_completions = set()
    print("Saving...")
    # for model in optimal_completions:
    #     one_min_sol_lst = set(a[0] for pred in model if pred == 'xreaction' for a in model[pred])
    #     all_completions.add(tuple(one_min_sol_lst))
    # all_completions = list(all_completions)
    with open(paths_map['results'], 'w') as f:
        f.write('execution time = {0}\n'.format(exec_time))
        f.write('objective value = {0}\n'.format(sol))
        f.write(f"reconstructable_targets = {reconstructable_targets}\n")
        f.write(f"unreconstructable = {never_producible}\n")
        f.write('minimal size = {0}\n'.format(optimum))
        f.write('One minimal completion of size {0}:\n\t{1}\n'.format(
            optimum, '\n\t'.join(optimum)))
        f.write('optimal completions:\n')
        # for completion in all_completions:
        #     f.write('\t{0}\n'.format('\n\t'.join(completion)))

def synechocystis():
    my_model = Model(model=read_sbml_model(prefix + 'tests/data/meneco/synechocystis/input/model.xml'), objective_function_id='e_Biomass__cytop')
    my_model.create_trnas_reactions()
    seeds, targets = [], []
    for metabolite in my_model.reactions.e_Biomass__cytop.reactants:
        targets.append((metabolite.id, metabolite.compartment))
    my_model.targets = targets
    my_model.identify_seeds()

    my_model.to_sbml('model', prefix + '/tests/data/meneco/synechocystis/input/', seeds=True, targets=True)

    paths_map = {'model': prefix + '/tests/data/meneco/synechocystis/input/model.xml',
                 'seeds': prefix + '/tests/data/meneco/synechocystis/input/model_seeds.xml',
                 'additional_seeds': prefix + '/tests/data/meneco/synechocystis/input/model_additional_seeds.xml',
                 'targets': prefix + '/tests/data/meneco/synechocystis/input/model_targets.xml',
                 'repairnet': prefix + '/tests/data/meneco/synechocystis/input/temporary_universal_model_transport.xml',
                 'results': prefix + '/tests/data/meneco/synechocystis/output/results.txt'}
    run_meneco(paths_map)


def lactis():
    # my_model = Model(model=read_sbml_model(prefix + 'tests/data/meneco/lactis/input/model.xml'), objective_function_id='e_Biomass__cytop')
    # my_model.create_trnas_reactions()
    # seeds, targets = [], []
    # for metabolite in my_model.reactions.e_Biomass__cytop.reactants:
    #     targets.append((metabolite.id, metabolite.compartment))
    # my_model.targets = targets
    # my_model.identify_seeds()
    # my_model.to_sbml('model', '../tests/data/meneco/lactis/input/', seeds=True, targets=True)

    paths_map = {'model': prefix + 'tests/data/meneco/lactis/input/model.xml',
                 'seeds': prefix + 'tests/data/meneco/lactis/input/model_seeds.xml',
                 'additional_seeds': prefix + 'tests/data/meneco/lactis/input/model_additional_seeds.xml',
                 'targets': prefix + 'tests/data/meneco/lactis/input/model_targets.xml',
                 'repairnet': prefix + 'tests/data/meneco/lactis/input/temporary_universal_model_transport.xml',
                 'results': prefix + 'tests/data/meneco/lactis/output/results.txt'}
    run_meneco(paths_map)


def cvulgaris():
    my_model = Model(model=read_sbml_model('../tests/data/meneco/cvulgaris/input/model.xml'), objective_function_id='e_Biomass__cytop')
    universal = read_sbml_model('../tests/data/meneco/cvulgaris/input/universal_model_compartmentalized_transport.xml')
    my_model.create_trnas_reactions()
    seeds, targets = [], []
    for metabolite in my_model.reactions.e_Biomass__cytop.reactants:
        targets.append((metabolite.id, metabolite.compartment))
    my_model.targets = targets #+ [(reactant.id, reactant.compartment) for reactant in my_model.reactions.e_Pigment__chlo.reactants]
    my_model.identify_seeds()
    my_model.add_seeds([(e.id, e.compartment) for e in universal.metabolites if
                        e.id in ("C00002__er",
                                 "C00002__chlo",
                                 "C00002__cytop",
                                 "C00031__cytop",
                                 "C03024__chlo",
                                 "C01092__mito",
                                 "C19845__mito"
                                 )])
    #
    my_model.to_sbml('model', '../tests/data/meneco/cvulgaris/input/', seeds=True, targets=True)

    paths_map = {'model': '../tests/data/meneco/cvulgaris/input/model.xml',
                 'seeds': '../tests/data/meneco/cvulgaris/input/model_seeds.xml',
                 'additional_seeds': '../tests/data/meneco/cvulgaris/input/model_additional_seeds.xml',
                 'targets': '../tests/data/meneco/cvulgaris/input/model_targets.xml',
                 'repairnet': '../tests/data/meneco/cvulgaris/input/universal_model_compartmentalized_transport.xml',
                 'results': '../tests/data/meneco/cvulgaris/output/results.txt'}
    run_meneco(paths_map)


def add_transport_to_compartimentalized():
    compartimentalized_model = read_sbml_model("../tests/data/meneco/lactis/input/universal_model_compartmentalized.xml")
    temporary_model = read_sbml_model(r"../tests/data/meneco/lactis/input/temporary_universal_model.xml")
    to_add = []
    for reaction in temporary_model.reactions:
        if reaction.id.startswith("T_"):
            to_add.append(reaction)
    compartimentalized_model.add_reactions(to_add)
    write_sbml_model(compartimentalized_model, r"../tests/data/meneco/lactis/input/temporary_universal_model_transport.xml")

def temp():
    # my_model = Model(model=read_sbml_model(r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\lactis\input\model_fixed.xml"), objective_function_id='e_Biomass__cytop')
    # seeds, targets = [], []
    # for metabolite in my_model.reactions.e_Biomass__in.reactants:
    #     targets.append((metabolite.id, metabolite.compartment))
    # my_model.targets = targets
    # my_model.identify_seeds()
    #
    # my_model.to_sbml('model', r"C:\Users\Bisbii\Desktop", seeds=True, targets=True)

    sys.path.append(r"C:\Users\Bisbii\PythonProjects\GSMMutils\src")
    from GSMMutils.model.COBRAmodel import MyModel
    from cobra.io import read_sbml_model
    import warnings
    import logging

    logging.getLogger("cobra").setLevel(logging.WARNING)
    warnings.filterwarnings("ignore")
    logging.getLogger("cobra").setLevel(logging.CRITICAL)
    logging.getLogger("meneco").setLevel(logging.CRITICAL)

    from cobra.io import write_sbml_model

    pre_model = MyModel(r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_fixed.xml", "e_Biomass__cytop")
    metabolite = pre_model.metabolites.C00028__cytop
    universal_model = read_sbml_model(r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\temporary_universal_model.xml")
    for reaction in metabolite.reactions:
        #model = MyModel(r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_fixed.xml", "e_Biomass__cytop")
        model = pre_model.copy()
        print(reaction.id)
        #temp = [r_id for r_id in model.metabolites.C00075__cytop.reactions if reaction]  #
        # print(set([r.id for r in model.metabolites.C00097__cytop.reactions]) - set([r.id for r in temp]))
        model.remove_reactions([reaction.id])
        #model.reactions.R02135__cytop.bounds = (-1000, 0)
        write_sbml_model(model, r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_fixed_v2.xml")


        paths_map = {
                    # 'model': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\lactis\input\model_fixed_v2.xml",
                     'model': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_fixed_v2.xml",
                     'seeds': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_seeds.xml",
                     'targets': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_targets - Copy.xml",
                     'repairnet': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\temporary_universal_model.xml"}

        draftnet = sbml.readSBMLnetwork(paths_map['model'], "draft")

        seeds = sbml.readSBMLseeds(paths_map['seeds'])

        targets = sbml.readSBMLtargets(paths_map['targets'])

        model = query.get_unproducible(draftnet, targets, seeds)
        unproducible = set([a[0] for pred in model if pred == 'unproducible_target' for a in model[pred]])
        if len(unproducible) > 0:
            print('{0} unproducible targets:\n\t{1}\n'.format(len(unproducible), '\n\t'.join(unproducible)))

        repairnet = sbml.readSBMLnetwork(paths_map['repairnet'], 'repairnet')

        combinet = draftnet
        combinet = TermSet(combinet.union(repairnet))

        model = query.get_unproducible(combinet, targets, seeds)
        never_producible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
        print('{0} unreconstructable targets:\n\t{1}\n'.format(len(never_producible), '\n\t'.join(never_producible)))

        reconstructable_targets = set(list(unproducible.difference(never_producible)))
        print('{0} reconstructable targets:\n\t{1}\n'.format(len(reconstructable_targets), '\n\t'.join(reconstructable_targets)))
        #

        #targets = TermSet(Atom('target("' + t + '")') for t in reconstructable_targets)
        optimum = extract_xreactions(query.get_minimal_completion_size(draftnet, repairnet, seeds, targets), False)

        print('minimal size =', optimum)

        print('One minimal completion of size {0}:\n\t{1}\n'.format(
            optimum, '\n\t'.join(optimum)))

        union = query.get_union_of_optimal_completions(draftnet, repairnet, seeds, targets, len(optimum))

        union = extract_xreactions(union, False)

        for r in optimum:
            in_model_reaction = universal_model.reactions.get_by_id(r.replace("R_",""))
            print(in_model_reaction.id, in_model_reaction.name)
            pre_model.add_reactions([in_model_reaction])
        write_sbml_model(pre_model, r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_fixed_v3.xml")

        optimal_completions = query.get_optimal_completions(draftnet, repairnet, seeds, targets, len(optimum))
        all_completions = set()
        for model in optimal_completions:
            one_min_sol_lst = set(a[0] for pred in model if pred == 'xreaction' for a in model[pred])
            all_completions.add(tuple(one_min_sol_lst))
        all_completions = list(all_completions)
        print("All completions: ", len(all_completions))
        print("All completions: ", all_completions)
        print(pre_model.slim_optimize())
        for completion in all_completions:
            temp_model = pre_model.copy()
            for r in completion:
                in_model_reaction = universal_model.reactions.get_by_id(r.replace("R_", ""))
                temp_model.add_reactions([in_model_reaction])
            print(temp_model.slim_optimize())


def temp2():
    sys.path.append(r"C:\Users\Bisbii\PythonProjects\GSMMutils\src")
    from GSMMutils.model.COBRAmodel import MyModel
    from cobra.io import read_sbml_model
    import warnings
    import logging


    mymodel = MyModel(r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_fixed.xml", "e_Biomass__cytop")
    # universal = read_sbml_model(r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\model_1\input\universal_model_compartmentalized.xml")
    mymodel.remove_reactions(["R03472__cytop"])
    # model.reactions.R02135__cytop.bounds = (-1000, 0)
    write_sbml_model(mymodel, r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_fixed_v2.xml")

    paths_map = {
        'model': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_fixed_v2.xml",
        'seeds': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_seeds.xml",
        'targets': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\model_targets.xml",
        'repairnet': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\synechocystis\input\temporary_universal_model.xml"}
    # paths_map = {
    #     'model': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\model_1\input\model_fixed_v3.xml",
    #     'seeds': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\model_1\input\model_seeds.xml",
    #     'targets': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\model_1\input\model_targets - Copy.xml",
    #     'repairnet': r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\model_1\input\universal_model_compartmentalized.xml"}

    draftnet = sbml.readSBMLnetwork(paths_map['model'], "draft")

    seeds = sbml.readSBMLseeds(paths_map['seeds'])

    targets = sbml.readSBMLtargets(paths_map['targets'])

    model = query.get_unproducible(draftnet, targets, seeds)
    unproducible = set([a[0] for pred in model if pred == 'unproducible_target' for a in model[pred]])
    if len(unproducible) > 0:
        print('{0} unproducible targets:\n\t{1}\n'.format(len(unproducible), '\n\t'.join(unproducible)))

    repairnet = sbml.readSBMLnetwork(paths_map['repairnet'], 'repairnet')

    combinet = draftnet
    combinet = TermSet(combinet.union(repairnet))

    model = query.get_unproducible(combinet, targets, seeds)
    never_producible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
    print('{0} unreconstructable targets:\n\t{1}\n'.format(len(never_producible), '\n\t'.join(never_producible)))

    reconstructable_targets = set(list(unproducible.difference(never_producible)))
    print('{0} reconstructable targets:\n\t{1}\n'.format(len(reconstructable_targets), '\n\t'.join(reconstructable_targets)))
    #

    # targets = TermSet(Atom('target("' + t + '")') for t in reconstructable_targets)
    optimum = extract_xreactions(query.get_minimal_completion_size(draftnet, repairnet, seeds, targets), False)

    print('minimal size =', optimum)

    print('One minimal completion of size {0}:\n\t{1}\n'.format(
        optimum, '\n\t'.join(optimum)))

    # for reaction in optimum:
    #     in_model_reaction = universal.reactions.get_by_id(reaction.replace("R_", ""))
    #     mymodel.add_reactions([in_model_reaction])
    # print(mymodel.slim_optimize())
    # write_sbml_model(mymodel, r"C:\Users\Bisbii\PythonProjects\gap_filling_dl\tests\data\model_1\input\model_fixed_v2.xml")


if __name__ == '__main__':
    # add_transport_to_compartimentalized()
    # lactis()
    # synechocystis()
    # cvulgaris()
    # start = time.time()
    # model = Model(read_sbml_model(prefix + r"/tests/data/meneco/cvulgaris/input/model.xml"), "e_Biomass__cytop")
    # print("Targets:\n")
    # targets = model.identify_targets()
    # print(targets)
    # print(len(targets))
    #temp()
    temp2()
    # print(time.time() - start)
    print("Done")
