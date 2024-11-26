from importlib import import_module
from inspect import isfunction

from jobflow import Flow, job
from jobflow.managers.local import run_locally
from rich.pretty import pprint
from universal_functions import get_input_dict, group_edges_dict, universal_reduce, get_workflow


def add_xy(x, y):
    return x + y, x, y


def add_xyz(x, y, z):
    return x + y + z


nodes_dict = {
    0: add_xy,
    1: add_xyz,
    2: 1,
    3: 2,
}


edges_lst = [
    {"target": 1, "targetHandle": "x", "source": 0, "sourceHandle": "x"},
    {"target": 1, "targetHandle": "y", "source": 0, "sourceHandle": "y"},
    {"target": 1, "targetHandle": "z", "source": 0, "sourceHandle": "z"},
    {"target": 0, "targetHandle": "x", "source": 2, "sourceHandle": None},
    {"target": 0, "targetHandle": "y", "source": 3, "sourceHandle": None},
]


nodes_updated_dict, edges_updated_lst = universal_reduce(
    nodes_dict=nodes_dict, edges_lst=edges_lst
)
nodes_updated_dict, edges_updated_lst


total_dict = group_edges_dict(edges_lst=edges_updated_lst)
input_dict = get_input_dict(nodes_dict=nodes_updated_dict)


task_lst = get_workflow(
    nodes_dict=nodes_updated_dict,
    input_dict=input_dict,
    total_dict=total_dict,
    backend='jobflow'
)

flow = Flow(task_lst)

result = run_locally(flow)
result
