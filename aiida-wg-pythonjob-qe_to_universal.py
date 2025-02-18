
# ! Specify Python Code manually to point to executable of my Python venv
# !


from pathlib import Path
import subprocess
from typing import Any
import shutil
import os
import ipdb
from ase import Atoms
from adis_tools.parsers import parse_pw
from ase.io import write

from aiida import load_profile, orm
from aiida_workgraph import WorkGraph, task
from atomistic_engine_unification.calculate_funcs_aiida_wg_pythonjob import (
    get_bulk_structure,
    calculate_qe,
    generate_structures,
)

load_profile()

pseudopotentials = {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"}
pseudo_path = Path(
    "/home/geiger_j/aiida_projects/adis/git-repos/compare-workflow-graphs/pseudos"
)
# pseudo_path = pathlib.Path.cwd() / "pseudos"
strain_lst = [0.9, 0.95, 1, 1.05, 1.10]

relax_input_dict = {
    "pseudopotentials": pseudopotentials,
    "kpts": (1, 1, 1),  # (3, 3, 3),
    # "kpts": (3, 3, 3),
    "calculation": "scf",  # "vc-relax",
    # "calculation": "vc-relax",
    "smearing": 0.02,
}

# @task.pythonjob(outputs=[{"name": "structure"}])
# def get_bulk_strcuture(): ...
get_bulk_structure_dec = task.pythonjob(outputs=[{"name": "structure"}])(
    get_bulk_structure
)


# @task.pythonjob(
#     outputs=[
#         {
#             "name": "scaled_atoms",
#             "identifier": "workgraph.namespace",
#             "metadata": {"dynamic": True},
#         }
#     ]
# )
# def generate_structures(): ...


generate_structures_dec = task.pythonjob(
    outputs=[
        {
            "name": "scaled_atoms",
            "identifier": "workgraph.namespace",
            "metadata": {"dynamic": True},
        }
    ]
)(generate_structures)


# @task.pythonjob(
#     outputs=[
#         {"name": "structure"},
#         {"name": "energy"},
#         {"name": "volume"},
#     ]
# )
# def calculate_qe(): ...

calculate_qe_dec = task.pythonjob(
    outputs=[
        {"name": "structure"},
        {"name": "energy"},
        {"name": "volume"},
    ]
)(calculate_qe)

#region
# @task.pythonjob(
#     outputs=[{"name": "result", "from": "context.result"}],
# )
# def all_scf(structures, input_dict):

# Outputs        PK    Type
# -------------  ----  ----------
# scaled_atoms
#     qe_0       3080  AtomsData
#     qe_1       3081  AtomsData
#     qe_2       3082  AtomsData
#     qe_3       3083  AtomsData
#     qe_4       3084  AtomsData
# remote_folder  3078  RemoteData
# retrieved      3079  FolderData

# Caller                 PK  Type
# -------------------  ----  ---------------
# generate_structures  3057  WorkGraph<test>

# for key, structure in structures.items():

#     relax_task = wg.add_task(
#         calculate_qe,
#         name=f"scf_{key}",
#         working_directory=orm.Str(key),
#         structure=structure,
#         input_dict=input_dict,
#     )

#     relax_task.set_context({f"result.{key}.structure": "structure"})
#     relax_task.set_context({f"result.{key}.energy": "energy"})
#     relax_task.set_context({f"result.{key}.volume": "volume"})

# return wg

# ipdb.set_trace()
#endregion

wg = WorkGraph("test")

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

# generate_structures_task = wg.add_task(
#     generate_structures_dec,
#     name="generate_structures",
#     structure=relax_task.outputs.structure,
#     strain_lst=strain_lst,
# )


wg.run()
