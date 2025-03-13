# ! Specify Python Code manually to point to executable of my Python venv


from pathlib import Path
import subprocess
from rich.pretty import pprint
from turtle import st
from typing import Any
import shutil
import os
import ipdb
from ase import Atoms
from adis_tools.parsers import parse_pw
from ase.io import write

from aiida import load_profile, orm
from aiida_workgraph import WorkGraph, task, map_
from atomistic_engine_unification.calculate_funcs_aiida_wg_pythonjob import (
    get_bulk_structure,
    calculate_qe,
    generate_structures,
    plot_energy_volume_curve,
    all_scf
)

load_profile()

pseudopotentials = {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"}
pseudo_path = Path(
    "/home/geiger_j/aiida_projects/adis/git-repos/compare-workflow-graphs/pseudos"
)
# pseudo_path = pathlib.Path.cwd() / "pseudos"
# strain_lst = [0.9, 0.95, 1, 1.05, 1.10]
strain_lst = [0.9, 1]

scf_input_dict = {
    "pseudopotentials": pseudopotentials,
    "kpts": (1, 1, 1),  # (3, 3, 3),
    # "kpts": (3, 3, 3),
    "calculation": "scf",  # "vc-relax",
    # "calculation": "vc-relax",
    "smearing": 0.02,
}


relax_input_dict = {
    "pseudopotentials": pseudopotentials,
    "kpts": (1, 1, 1),  # (3, 3, 3),
    # "kpts": (3, 3, 3),
    "calculation": "scf",  # "vc-relax",
    # "calculation": "vc-relax",
    "smearing": 0.02,
}

get_bulk_structure_dec = task.pythonjob(outputs=[{"name": "structure"}])(
    get_bulk_structure
)

generate_structures_dec = task.pythonjob(
    outputs=[
        {
            "name": "scaled_atoms",
            "identifier": "workgraph.namespace",
            "metadata": {"dynamic": True},
        }
    ]
)(generate_structures)

calculate_qe_dec = task.pythonjob(
    outputs=[
        {"name": "structure"},
        {"name": "energy"},
        {"name": "volume"},
    ]
)(calculate_qe)

all_scf_dec = task.pythonjob(
    outputs=[
        {
            "name": "qe_results",
            # "identifier": "workgraph.namespace",
            # "metadata": {"dynamic": True},
        }
    ]
)(all_scf)

plot_energy_volume_curve_dec = task.pythonjob()(plot_energy_volume_curve)


wg = WorkGraph("wg-pythonjob")

get_bulk_structure_task = wg.add_task(
    get_bulk_structure_dec,
    name="get_bulk_structure",
    element="Al",
    a=4.05,
    cubic=True,
)

relax_task = wg.add_task(
    calculate_qe_dec,
    name="relax",
    structure=get_bulk_structure_task.outputs.structure,
    input_dict=relax_input_dict,
    working_directory="relax",
)

generate_structures_task = wg.add_task(
    generate_structures_dec,
    name="generate_structures",
    structure=relax_task.outputs.structure,
    strain_lst=strain_lst,
)

map_(wg.tasks.generate_structures.outputs.scaled_atoms)(
    wg.add_task(
        calculate_qe_dec,
        name='scf',
        structure='{{map_input}}',
        working_directory='scf',  # scf_key
        input_dict=scf_input_dict,
    ),
    # wg.tasks.scf.set(scf_input_dict),
)

# all_scf_task = wg.add_task(
#     all_scf_dec,
#     name='all_scf',
#     structures=generate_structures_task.outputs.scaled_atoms,
#     input_dict=scf_input_dict
# )

wg.to_html('wg-pythonjob_map.html')
wg_dict = wg.to_dict()
# wg.run()

# plot_energy_volume_curve_task = wg.add_task(
#     plot_energy_volume_curve_dec,
#     name='plot_energy_volume_curve',
#     qe_results=all_scf_task.outputs.qe_results,
# )

# import ipdb; ipdb.set_trace()

# raise SystemExit()

# wg.run()


# ! This seems eerily similar to the `NodeLink` class from node-graph
# def get_edges_list(wg_dict):

#     edges_label_lst = []
#     for link_dict in wg_dict["links"]:
#         if link_dict["from_socket"] == "result":
#             edges_label_lst.append(
#                 {
#                     "target": link_dict["to_node"],
#                     "targetHandle": link_dict["to_socket"],
#                     "source": link_dict["from_node"],
#                     "sourceHandle": None,
#                 }
#             )
#         else:
#             edges_label_lst.append(
#                 {
#                     "target": link_dict["to_node"],
#                     "targetHandle": link_dict["to_socket"],
#                     "source": link_dict["from_node"],
#                     "sourceHandle": link_dict["from_socket"],
#                 }
#             )

#     return edges_label_lst

# edges_label_lst = get_edges_list(wg_dict=wg_dict)

# # print("EDGES_LABEL_LIST")
# # pprint(edges_label_lst)


# kwargs_dict, function_dict = {}, {}

# from node_graph.executor import NodeExecutor

# # for task_name, task_dict in wg_dict["tasks"].items():

# #     # Extract relevant input variable names (excluding metadata and _wait)
# #     input_variables = []
# #     for param in task_dict["inputs"]:
# #         if not param.startswith("metadata") and not param.startswith("_wait") and not param.startswith('monitors'):
# #             input_variables.append(param)

# #     # Prepare input keyword arguments
# #     input_kwargs = {}

# #     for param in input_variables:
# #         print(f"PARAM: {param}")
# #         try:
# #             property_value = task_dict["inputs"][param]["property"]["value"]
# #         except:
# #             # ipdb.set_trace()
# #             # raise
# #             pass

# #         if isinstance(property_value, dict):
# #             try:
# #                 input_kwargs[param] = property_value.value
# #             except:
# #                 input_kwargs[param] = property_value
# #                 # ipdb.set_trace()
# #         else:
# #             input_kwargs[param] = property_value

# #     # function_dict[task_name] = loads(task_dict['executor']['callable']).process_class._func
# #     kwargs_dict[task_name] = input_kwargs
# #     try:
# #         executor = NodeExecutor(**task_dict['executor']).executor
# #         function_dict[task_name] = executor
# #     except:

# #         ipdb.set_trace()

# #     # function_dict[task_name] = task_dict["executor"]["callable"]


# # print("KWARGS_DICT")
# # pprint(kwargs_dict)
# # print("FUNCTION_DICT")
# # pprint(function_dict)

# def write_workflow_json(wg):
#     wgdata = wg.to_dict()
#     data = {"nodes": {}, "edges": []}
#     node_name_mapping = {}
#     i = 0
#     for name, node in wgdata["tasks"].items():
#         ipdb.set_trace()
#         node_name_mapping[name] = i
#         callable_name = node["executor"]["callable_name"]
#         data["nodes"][i] = callable_name
#         if callable_name == "pickle_node":
#             data["nodes"][i] = node["inputs"]["value"]["property"]["value"].value
#         i += 1

#     for link in wgdata["links"]:
#         if wgdata["tasks"][link["from_node"]]["executor"]["callable_name"] == "pickle_node":
#             link["from_socket"] = None
#         link["from_node"] = node_name_mapping[link["from_node"]]
#         link["to_node"] = node_name_mapping[link["to_node"]]
#         data["edges"].append(link)

#     return data

# wf_data = write_workflow_json(wg=wg)
# pprint(wf_data)