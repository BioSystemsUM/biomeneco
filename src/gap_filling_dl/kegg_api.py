import json

from KEGG_parser.parsers import parse_pathway
from requests import get
from Bio.KEGG.REST import kegg_get


def search_pathway_map_id(pathway_name):
    url = f'https://rest.kegg.jp/find/pathway/{pathway_name}'
    response = get(url)
    for line in response.text.splitlines():
        if len(line.split('\t')) == 1:
            continue
        pathway_id, pathway_name = line.split('\t')
        if pathway_name == pathway_name:
            return pathway_id.strip("path:")


def get_kegg_pathways(pathway_id=None):
    """Get KEGG pathways from KEGG API"""
    if pathway_id:
        try:
            result = kegg_get(pathway_id)
        except Exception as e:
            print(pathway_id)
            print(e)
            return None
        if result:
            result = result.read()
        result = parse_pathway(result.replace("///", ""))
        return result


def create_related_pathways_map(universal_model, folder_path):
    res = {}
    for pathway in universal_model.groups:
        related_pathways = get_related_pathways(pathway.name, res)
        res[pathway.name] = related_pathways
    json.dump(res, open(f'{folder_path}/../related_pathways_map.json', 'w'))
    return res


def get_related_pathways(pathway_name, related_pathways_map):
    if pathway_name in related_pathways_map:
        return related_pathways_map[pathway_name]
    pathways = []
    try:
        pathway_id = search_pathway_map_id(pathway_name)
        if pathway_id is not None:
            kegg_result = get_kegg_pathways(pathway_id)
            for pathway in kegg_result['REL_PATHWAY']:
                pathways.append(pathway[1])
    except Exception as e:
        print(e)
    return pathways
