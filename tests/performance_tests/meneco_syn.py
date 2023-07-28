from cobra.io import read_sbml_model
from clyngor.as_pyasp import TermSet,Atom
from urllib.request import urlopen
from meneco.meneco import query, utils, sbml

import logging

logging.getLogger("cobra").setLevel(logging.WARNING)

draftnet = sbml.readSBMLnetwork('../data/draft_synechocystis/input/model.xml', "draftnet")


seeds = sbml.readSBMLseeds('../data/draft_synechocystis/draft_synechocystis_seeds.xml')


targets = sbml.readSBMLtargets('../data/draft_synechocystis/draft_synechocystis_targets.xml')

model = query.get_unproducible(draftnet, targets, seeds)
unproducible = set([a[0] for pred in model if pred == 'unproducible_target' for a in model[pred]])
print('{0} unproducible targets:\n\t{1}\n'.format(len(unproducible), '\n\t'.join(unproducible)))


repairnet = sbml.readSBMLnetwork('../data/draft_synechocystis/temporary_universal_model.xml', 'repairnet')
# repairnet = sbml.readSBMLnetwork('../data/draft_synechocystis/universal_model_compartmentalized.xml', 'repairnet')

combinet = draftnet
combinet = TermSet(combinet.union(repairnet))

model = query.get_unproducible(combinet, targets, seeds)
never_producible = set(a[0] for pred in model if pred == 'unproducible_target' for a in model[pred])
print('{0} unreconstructable targets:\n\t{1}\n'.format(len(never_producible), '\n\t'.join(never_producible)))

reconstructable_targets = set(list(unproducible.difference(never_producible)))
print('{0} reconstructable targets:\n\t{1}\n'.format(len(reconstructable_targets), '\n\t'.join(reconstructable_targets)))
#
# essential_reactions = set()
# for t in reconstructable_targets:
#       single_target = TermSet()
#       single_target.add(Atom('target("' + t + '")'))
#       print('\nComputing essential reactions for', t,'... ', end=' ')
#       model = query.get_intersection_of_completions(draftnet, repairnet, seeds, single_target)
#       print(' done.')
#       essentials_lst = set(a[0] for pred in model if pred == 'xreaction' for a in model[pred])
#       print('{0} essential reactions for target {1}:\n\t{2}'.format(len(essentials_lst), t, '\n\t'.join(essentials_lst)))
#       essential_reactions = essential_reactions.union(essentials_lst)
# print('Overall {0} essential reactions found:\n\t{1}\n'.format(len(essential_reactions), '\n\t'.join(essential_reactions)))

targets = TermSet(Atom('target("' + t+'")') for t in reconstructable_targets)
model = query.get_minimal_completion_size(draftnet, repairnet, seeds, targets)
one_min_sol_lst = set(a[0] for pred in model if pred == 'xreaction' for a in model[pred])
optimum = len(one_min_sol_lst)

print('minimal size =',optimum)

print('One minimal completion of size {0}:\n\t{1}\n'.format(
            optimum, '\n\t'.join(one_min_sol_lst)))

with open('results.txt', 'w') as f:
    f.write('minimal size = {0}\n'.format(optimum))
    f.write('One minimal completion of size {0}:\n\t{1}\n'.format(
            optimum, '\n\t'.join(one_min_sol_lst)))
