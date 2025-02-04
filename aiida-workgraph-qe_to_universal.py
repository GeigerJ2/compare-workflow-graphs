import os
import pathlib
import shutil
import subprocess

import matplotlib.pyplot as plt
from adis_tools.parsers import parse_pw
from aiida import load_profile, orm
from aiida_workgraph import WorkGraph, task
from ase.atoms import Atoms
from ase.build import bulk
from ase.io import write

load_profile()

pseudopotentials = {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"}
# pseudo_path = pathlib.Path(
#     "/home/geiger_j/aiida_projects/adis/git-repos/compare-workflow-graphs/pseudos"
# )
pseudo_path = pathlib.Path.cwd() / "pseudos"
strain_lst = orm.List([0.9, 0.95, 1, 1.05, 1.10])

relax_input_dict = orm.Dict(
    {
        "pseudopotentials": pseudopotentials,
        # "kpts": (1, 1, 1),  # (3, 3, 3),
        "kpts": (3, 3, 3),
        # "calculation": "scf",  # "vc-relax",
        "calculation": "vc-relax",
        "smearing": 0.02,
    }
)

scf_input_dict = orm.Dict(
    {
        "pseudopotentials": pseudopotentials,
        "kpts": (3, 3, 3),
        "calculation": "scf",  # "vc-relax",
        "smearing": 0.02,
    }
)


@task.calcfunction(outputs=[{"name": "structure"}])
def get_bulk_structure(element: orm.Str, a: orm.Float, cubic: orm.Bool):
    atoms = bulk(name=element.value, a=a.value, cubic=cubic.value)
    return {"structure": orm.StructureData(ase=atoms)}


# ! Cannot use integers as link labels, as they are not valid python identifiers
@task.calcfunction(outputs=[{"name": "structures"}])
def generate_structures(structure: orm.StructureData, strain_lst: orm.List):
    structure_lst = []
    for strain in strain_lst:
        atoms = structure.get_ase()
        structure_strain = atoms.copy()
        structure_strain.set_cell(
            structure_strain.cell * strain ** (1 / 3), scale_atoms=True
        )
        structure_lst.append(structure_strain)

    return_dict = {
        f"qe_{str(i)}": orm.StructureData(ase=atoms)
        for i, atoms in enumerate(structure_lst)
    }
    return {"structures": return_dict}


@task.calcfunction(
    outputs=[{"name": "structure"}, {"name": "energy"}, {"name": "volume"}]
)
def collect_output(working_directory: orm.Str):
    output = parse_pw(os.path.join(working_directory.value, "pwscf.xml"))
    output_orm_dict = {
        "structure": orm.StructureData(ase=output["ase_structure"]),
        "energy": orm.Float(output["energy"]),
        "volume": orm.Float(output["ase_structure"].get_volume()),
    }
    return output_orm_dict


@task.calcfunction()
def write_input(
    input_dict: orm.Dict, structure: orm.StructureData, working_directory: orm.Str
):
    filename = os.path.join(working_directory.value, "input.pwi")
    os.makedirs(working_directory.value, exist_ok=True)
    atoms = structure.get_ase()
    python_dict = input_dict.get_dict()
    pseudopotentials = python_dict["pseudopotentials"]
    shutil.copy(src=pseudo_path / pseudopotentials["Al"], dst=working_directory.value)
    write(
        filename=filename,
        images=atoms,
        Crystal=True,
        kpts=python_dict["kpts"],
        input_data={
            "calculation": python_dict["calculation"],
            "occupations": "smearing",
            "degauss": python_dict["smearing"],
            "pseudo_dir": "./",
        },
        pseudopotentials=pseudopotentials,
        tstress=True,
        tprnfor=True,
    )


@task.calcfunction()
def run_qe(working_directory: orm.Str):
    _ = subprocess.check_output(
        "mpirun -np 1 pw.x -in input.pwi > output.pwo",
        cwd=working_directory.value,
        shell=True,
    )


# ! Or, maybe, just make the whole process of writing the input and running the QE process as one calcfunction
@task.graph_builder(
    outputs=[
        {"name": "structure", "from": "collect_output.structure"},
        {"name": "energy", "from": "collect_output.energy"},
        {"name": "volume", "from": "collect_output.volume"},
    ]
)
def calculate_qe(
    working_directory: orm.Str, input_dict: orm.Dict, structure: orm.StructureData
):
    wg = WorkGraph("calculate_qe")
    write_input_task = wg.add_task(
        write_input,
        name="write_input",
        input_dict=input_dict,
        structure=structure,
        working_directory=working_directory,
    )
    run_qe_task = wg.add_task(
        run_qe,
        name="run_qe",
        working_directory=working_directory,
    )
    run_qe_task.waiting_on.add(write_input_task)

    collect_output_task = wg.add_task(
        collect_output, name="collect_output", working_directory=working_directory
    )
    collect_output_task.waiting_on.add(run_qe_task)

    return wg


@task.graph_builder(
    outputs=[{"name": "result", "from": "context.result"}],
)
def all_scf(structures, input_dict):
    wg = WorkGraph("all_scf")
    for key, structure in structures.items():
        relax_task = wg.add_task(
            calculate_qe,
            name=f"scf_{key}",
            working_directory=orm.Str(key),
            structure=structure,
            input_dict=input_dict,
        )

        relax_task.set_context({f"result.{key}.structure": "structure"})
        relax_task.set_context({f"result.{key}.energy": "energy"})
        relax_task.set_context({f"result.{key}.volume": "volume"})

    return wg


@task.calcfunction()
def plot_energy_volume_curve(**datas):
    volume_lst, energy_lst = [], []
    for key, data in datas.items():
        volume_lst.append(data["volume"].value)
        energy_lst.append(data["energy"].value)

    plt.plot(volume_lst, energy_lst)
    plt.xlabel("Volume")
    plt.ylabel("Energy")
    plt.savefig("evcurve.png")


wg_eos = WorkGraph("eos")

get_bulk_structure_task = wg_eos.add_task(
    get_bulk_structure,
    name="get_bulk_structure",
    element=orm.Str("Al"),
    a=orm.Float(4.05),
    cubic=orm.Bool(True),
)

relax_task = wg_eos.add_task(
    calculate_qe,
    name="calculate_qe",
    working_directory=orm.Str("relax"),
    structure=get_bulk_structure_task.outputs.structure,
    input_dict=relax_input_dict,
)

scale_task = wg_eos.add_task(
    generate_structures,
    name="generate_structures",
    structure=relax_task.outputs.structure,
    strain_lst=strain_lst,
)

all_scf_task = wg_eos.add_task(
    all_scf,
    name="all_scf",
    structures=scale_task.outputs.structures,
    input_dict=scf_input_dict,
)

plot_task = wg_eos.add_task(plot_energy_volume_curve, datas=all_scf_task.outputs.result)

wg_eos.run()
