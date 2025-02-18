
from pyiron_base import Project

from from_universal_funcs_pyiron import (
    get_delayed_object_dict,
    get_dict,
    # get_kwargs,
    get_list,
    # get_source,
    get_source_handles,
    group_edges,
    resort_total_lst,
)

from to_universal_funcs_pyiron import (
    calculate_qe,
    get_bulk_structure,
    generate_structures,
    plot_energy_volume_curve,
)


pseudopotentials = {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"}


nodes_new_dict = {  # from jobflow
    0: get_bulk_structure,
    1: calculate_qe,
    2: generate_structures,
    3: calculate_qe,
    4: calculate_qe,
    5: calculate_qe,
    6: calculate_qe,
    7: calculate_qe,
    8: plot_energy_volume_curve,
    9: "Al",
    10: 4.05,
    11: True,
    12: "mini",
    13: get_dict,
    14: {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
    15: [3, 3, 3],
    16: "vc-relax",
    17: 0.02,
    18: [0.9, 0.9500000000000001, 1.0, 1.05, 1.1],
    19: "strain_0",
    20: get_dict,
    21: "scf",
    22: "strain_1",
    23: get_dict,
    24: "strain_2",
    25: get_dict,
    26: "strain_3",
    27: get_dict,
    28: "strain_4",
    29: get_dict,
    30: get_list,
    31: get_list,
}


edges_new_lst = [
    {"target": 0, "targetHandle": "name", "source": 9, "sourceHandle": None},
    {"target": 0, "targetHandle": "a", "source": 10, "sourceHandle": None},
    {"target": 0, "targetHandle": "cubic", "source": 11, "sourceHandle": None},
    {
        "target": 1,
        "targetHandle": "working_directory",
        "source": 12,
        "sourceHandle": None,
    },
    {"target": 13, "targetHandle": "structure", "source": 0, "sourceHandle": None},
    {
        "target": 13,
        "targetHandle": "pseudopotentials",
        "source": 14,
        "sourceHandle": None,
    },
    {"target": 13, "targetHandle": "kpts", "source": 15, "sourceHandle": None},
    {"target": 13, "targetHandle": "calculation", "source": 16, "sourceHandle": None},
    {"target": 13, "targetHandle": "smearing", "source": 17, "sourceHandle": None},
    {"target": 1, "targetHandle": "input_dict", "source": 13, "sourceHandle": None},
    {
        "target": 2,
        "targetHandle": "structure",
        "source": 1,
        "sourceHandle": "structure",
    },
    {"target": 2, "targetHandle": "strain_lst", "source": 18, "sourceHandle": None},
    {
        "target": 3,
        "targetHandle": "working_directory",
        "source": 19,
        "sourceHandle": None,
    },
    {"target": 20, "targetHandle": "structure", "source": 2, "sourceHandle": "0"},
    {
        "target": 20,
        "targetHandle": "pseudopotentials",
        "source": 14,
        "sourceHandle": None,
    },
    {"target": 20, "targetHandle": "kpts", "source": 15, "sourceHandle": None},
    {"target": 20, "targetHandle": "calculation", "source": 21, "sourceHandle": None},
    {"target": 20, "targetHandle": "smearing", "source": 17, "sourceHandle": None},
    {"target": 3, "targetHandle": "input_dict", "source": 20, "sourceHandle": None},
    {
        "target": 4,
        "targetHandle": "working_directory",
        "source": 22,
        "sourceHandle": None,
    },
    {"target": 23, "targetHandle": "structure", "source": 2, "sourceHandle": "1"},
    {
        "target": 23,
        "targetHandle": "pseudopotentials",
        "source": 14,
        "sourceHandle": None,
    },
    {"target": 23, "targetHandle": "kpts", "source": 15, "sourceHandle": None},
    {"target": 23, "targetHandle": "calculation", "source": 21, "sourceHandle": None},
    {"target": 23, "targetHandle": "smearing", "source": 17, "sourceHandle": None},
    {"target": 4, "targetHandle": "input_dict", "source": 23, "sourceHandle": None},
    {
        "target": 5,
        "targetHandle": "working_directory",
        "source": 24,
        "sourceHandle": None,
    },
    {"target": 25, "targetHandle": "structure", "source": 2, "sourceHandle": "2"},
    {
        "target": 25,
        "targetHandle": "pseudopotentials",
        "source": 14,
        "sourceHandle": None,
    },
    {"target": 25, "targetHandle": "kpts", "source": 15, "sourceHandle": None},
    {"target": 25, "targetHandle": "calculation", "source": 21, "sourceHandle": None},
    {"target": 25, "targetHandle": "smearing", "source": 17, "sourceHandle": None},
    {"target": 5, "targetHandle": "input_dict", "source": 25, "sourceHandle": None},
    {
        "target": 6,
        "targetHandle": "working_directory",
        "source": 26,
        "sourceHandle": None,
    },
    {"target": 27, "targetHandle": "structure", "source": 2, "sourceHandle": "3"},
    {
        "target": 27,
        "targetHandle": "pseudopotentials",
        "source": 14,
        "sourceHandle": None,
    },
    {"target": 27, "targetHandle": "kpts", "source": 15, "sourceHandle": None},
    {"target": 27, "targetHandle": "calculation", "source": 21, "sourceHandle": None},
    {"target": 27, "targetHandle": "smearing", "source": 17, "sourceHandle": None},
    {"target": 6, "targetHandle": "input_dict", "source": 27, "sourceHandle": None},
    {
        "target": 7,
        "targetHandle": "working_directory",
        "source": 28,
        "sourceHandle": None,
    },
    {"target": 29, "targetHandle": "structure", "source": 2, "sourceHandle": "4"},
    {
        "target": 29,
        "targetHandle": "pseudopotentials",
        "source": 14,
        "sourceHandle": None,
    },
    {"target": 29, "targetHandle": "kpts", "source": 15, "sourceHandle": None},
    {"target": 29, "targetHandle": "calculation", "source": 21, "sourceHandle": None},
    {"target": 29, "targetHandle": "smearing", "source": 17, "sourceHandle": None},
    {"target": 7, "targetHandle": "input_dict", "source": 29, "sourceHandle": None},
    {"target": 30, "targetHandle": "0", "source": 3, "sourceHandle": "volume"},
    {"target": 30, "targetHandle": "1", "source": 4, "sourceHandle": "volume"},
    {"target": 30, "targetHandle": "2", "source": 5, "sourceHandle": "volume"},
    {"target": 30, "targetHandle": "3", "source": 6, "sourceHandle": "volume"},
    {"target": 30, "targetHandle": "4", "source": 7, "sourceHandle": "volume"},
    {"target": 8, "targetHandle": "volume_lst", "source": 30, "sourceHandle": None},
    {"target": 31, "targetHandle": "0", "source": 3, "sourceHandle": "energy"},
    {"target": 31, "targetHandle": "1", "source": 4, "sourceHandle": "energy"},
    {"target": 31, "targetHandle": "2", "source": 5, "sourceHandle": "energy"},
    {"target": 31, "targetHandle": "3", "source": 6, "sourceHandle": "energy"},
    {"target": 31, "targetHandle": "4", "source": 7, "sourceHandle": "energy"},
    {"target": 8, "targetHandle": "energy_lst", "source": 31, "sourceHandle": None},
]


pr = Project("test")
pr.remove_jobs(recursive=True, silently=True)


total_lst = group_edges(edges_new_lst)


total_new_lst = resort_total_lst(total_lst=total_lst, nodes_dict=nodes_new_dict)


source_handle_dict = get_source_handles(edges_new_lst)


delayed_object_dict = get_delayed_object_dict(
    total_lst=total_new_lst,
    nodes_dict=nodes_new_dict,
    source_handle_dict=source_handle_dict,
    pyiron_project=pr,
)


delayed_object_dict[list(delayed_object_dict.keys())[-1]].draw()


delayed_object_dict[list(delayed_object_dict.keys())[-1]].pull()


pr.job_table()
