#!/usr/bin/env python
# coding: utf-8


from pyiron_base import Project, job
from pyiron_base.project.delayed import DelayedObject
from universal_functions import get_kwargs, group_edges_list, get_source_handles


from inspect import isfunction


def get_source(nodes_dict, delayed_object_dict, source, sourceHandle):
    if source in delayed_object_dict.keys():
        return (
            delayed_object_dict[source].__getattr__("output").__getattr__(sourceHandle)
        )
    else:
        return nodes_dict[source]


def get_delayed_object_dict(total_lst, nodes_dict, source_handle_dict, pyiron_project):
    delayed_object_dict = {}
    for item in total_lst:
        key, input_dict = item
        kwargs = {
            k: get_source(
                nodes_dict=nodes_dict,
                delayed_object_dict=delayed_object_dict,
                source=v["source"],
                sourceHandle=v["sourceHandle"],
            )
            for k, v in input_dict.items()
        }
        delayed_object_dict[key] = job(
            funct=nodes_dict[key],
            output_key_lst=source_handle_dict.get(key, []),
        )(**kwargs, pyiron_project=pyiron_project)
    return delayed_object_dict

def add_x_and_y(x, y):
    z = x + y
    return {"x": x, "y": y, "z": z}


def add_x_and_y_and_z(x, y, z):
    w = x + y + z
    return w


nodes_dict = {
    0: add_x_and_y_and_z,
    1: add_x_and_y,
    2: 1,
    3: 2,
}


edges_lst = [
    {"target": 0, "targetHandle": "x", "source": 1, "sourceHandle": "x"},
    {"target": 1, "targetHandle": "x", "source": 2, "sourceHandle": None},
    {"target": 1, "targetHandle": "y", "source": 3, "sourceHandle": None},
    {"target": 0, "targetHandle": "y", "source": 1, "sourceHandle": "y"},
    {"target": 0, "targetHandle": "z", "source": 1, "sourceHandle": "z"},
]


nodes_dict, edges_lst


pr = Project("test")


pr.remove_jobs(recursive=True, silently=True)


delayed_object_dict = get_delayed_object_dict(
    total_lst=group_edges_list(edges_lst),
    nodes_dict=nodes_dict,
    source_handle_dict=get_source_handles(edges_lst),
    pyiron_project=pr,
)
delayed_object_dict


delayed_object_dict[list(delayed_object_dict.keys())[-1]].pull()


pr.job_table()
