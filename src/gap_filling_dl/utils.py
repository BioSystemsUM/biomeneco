import json
import os

import cobra
from cobra.io import write_sbml_model, read_sbml_model

from biomeneco.model import Model
from kegg_api import get_related_pathways, create_related_pathways_map



