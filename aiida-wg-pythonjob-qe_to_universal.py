# ! Specify Python Code manually to point to executable of my Python venv


from pathlib import Path
import subprocess
from turtle import st
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
    plot_energy_volume_curve,
    all_scf
)

load_profile()

pseudopotentials = {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"}
pseudo_path = Path(
    "/home/geiger_j/aiida_projects/adis/git-repos/compare-workflow-graphs/pseudos"
)
# pseudo_path = pathlib.Path.cwd() / "pseudos"
strain_lst = [0.9, 0.95, 1, 1.05, 1.10]
# strain_lst = [0.9, 1]

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
    # "calculation": "scf",  # "vc-relax",
    "calculation": "vc-relax",
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

# TODO: Try also just normal code here, appending to WG within the for-loop
# TODO: Also try to build an all_scf WG and pass it as a task
# TODO: `return_atoms` optional
# TODO: Convert Atoms into pure dict, also replacing np.bool

all_scf_dec = task.pythonjob(
    outputs=[
        {
            "name": "qe_results",
            # "identifier": "workgraph.namespace",
            # "metadata": {"dynamic": True},
        }
    ]
)(all_scf)


all_scf_task = wg.add_task(
    all_scf_dec,
    name='all_scf',
    structures=generate_structures_task.outputs.scaled_atoms,
    input_dict=scf_input_dict
)

plot_energy_volume_curve_task = wg.add_task(
    plot_energy_volume_curve_dec,
    name='plot_energy_volume_curve',
    qe_results=all_scf_task.outputs.qe_results,
)

# d = wg.to_dict()

# import ipdb; ipdb.set_trace()

wg.run()
