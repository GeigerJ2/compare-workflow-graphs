import ipdb
import numpy as np
from pyiron_base import Project, job

# from .calculate_funcs_pyironcalculate_funcs_pyiron import (
#     calculate_qe,
#     generate_structures,
#     get_bulk_structure,
#     plot_energy_volume_curve,
# )
import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from adis_tools.parsers import parse_pw
from ase.atoms import Atoms
from ase.build import bulk
from ase.io import write
from pyiron_base import Project, job
from pyiron_base.project.delayed import DelayedObject, draw


def write_input(input_dict, working_directory="."):
    filename = os.path.join(working_directory, "input.pwi")
    os.makedirs(working_directory, exist_ok=True)
    write(
        filename=filename,
        images=Atoms(**input_dict["structure"]),
        Crystal=True,
        kpts=input_dict["kpts"],
        input_data={
            "calculation": input_dict["calculation"],
            "occupations": "smearing",
            "degauss": input_dict["smearing"],
        },
        pseudopotentials=input_dict["pseudopotentials"],
        tstress=True,
        tprnfor=True,
    )


def collect_output(working_directory="."):
    output = parse_pw(os.path.join(working_directory, "pwscf.xml"))
    return {
        "structure": output["ase_structure"].todict(),
        "energy": output["energy"],
        "volume": output["ase_structure"].get_volume(),
    }


def calculate_qe(working_directory, input_dict):
    write_input(
        input_dict=input_dict,
        working_directory=working_directory,
    )
    subprocess.check_output(
        "mpirun -np 1 pw.x -in input.pwi > output.pwo",
        cwd=working_directory,
        shell=True,
    )
    return collect_output(working_directory=working_directory)


def get_bulk_structure(name, a, cubic):
    return bulk(
        name=name,
        a=a,
        cubic=cubic,
    ).todict()


def generate_structures(structure, strain_lst):
    structure_lst = []
    for strain in strain_lst:
        structure_strain = Atoms(**structure)
        structure_strain.set_cell(
            structure_strain.cell * strain ** (1 / 3), scale_atoms=True
        )
        structure_lst.append(structure_strain)
    return {str(i): s.todict() for i, s in enumerate(structure_lst)}


def plot_energy_volume_curve(volume_lst, energy_lst):
    plt.plot(volume_lst, energy_lst)
    plt.xlabel("Volume")
    plt.ylabel("Energy")
    plt.savefig("evcurve.png")


from to_universal_funcs_pyiron import (
    get_connection_dict,
    get_nodes,
    get_unique_objects,
    remove_server_obj,
)

pseudopotentials = {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"}


calculate_qe = job(output_key_lst=["energy", "volume", "structure"])(calculate_qe)
generate_structures = job()(generate_structures)
plot_energy_volume_curve = job()(plot_energy_volume_curve)


pr = Project(".")
pr.remove_jobs(recursive=True, silently=True)


structure = get_bulk_structure(
    name="Al",
    a=4.05,
    cubic=True,
    # pyiron_project=pr,
)


calc_mini = calculate_qe(
    working_directory="mini",
    input_dict={
        "structure": structure,
        "pseudopotentials": pseudopotentials,
        "kpts": (3, 3, 3),
        "calculation": "vc-relax",
        "smearing": 0.02,
    },
    # pyiron_project=pr,
)


number_of_strains = 5

structure_lst = generate_structures(  # the generate_structures() function is not available in the workflow graph
    structure=calc_mini.output.structure,
    strain_lst=np.linspace(0.9, 1.1, number_of_strains),
    # pyiron_project=pr,
    list_length=number_of_strains,
)

import ipdb; ipdb.set_trace()

job_strain_lst = []
for i, structure_strain in enumerate(structure_lst):
    calc_strain = calculate_qe(
        working_directory="strain_" + str(i),
        input_dict={
            "structure": structure_strain,
            "pseudopotentials": pseudopotentials,
            "kpts": (3, 3, 3),
            "calculation": "scf",
            "smearing": 0.02,
        },
        # pyiron_project=pr,
    )
    job_strain_lst.append(calc_strain)


plot = plot_energy_volume_curve(
    volume_lst=[job.output.volume for job in job_strain_lst],
    energy_lst=[job.output.energy for job in job_strain_lst],
    # pyiron_project=pr,
)


# This concludes the first version of the simulation workflow, in the following the submission to HPC resources, the different options for data storage and the publication of the workflow are briefly discussed.


plot.pull()


plot.draw()


nodes_dict, edges_lst = plot.get_graph()


nodes_dict, edges_lst = remove_server_obj(nodes_dict=nodes_dict, edges_lst=edges_lst)


delayed_object_updated_dict, match_dict = get_unique_objects(
    nodes_dict=nodes_dict, edges_lst=edges_lst
)


connection_dict, lookup_dict = get_connection_dict(
    delayed_object_updated_dict=delayed_object_updated_dict, match_dict=match_dict
)


nodes_new_dict = get_nodes(
    connection_dict=connection_dict,
    delayed_object_updated_dict=delayed_object_updated_dict,
)
nodes_new_dict
