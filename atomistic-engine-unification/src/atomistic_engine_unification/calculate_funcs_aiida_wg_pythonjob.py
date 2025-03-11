from typing import TYPE_CHECKING
import ipdb
from ase import Atoms  # This is only in type hints, so this import should not be necessary

if TYPE_CHECKING:
    from pathlib import Path
    from ase import Atoms


__all__ = (
    "get_bulk_structure",
    "generate_structures",
    "calculate_qe",
    "all_scf",
    "plot_energy_volume_curve",
)


def get_bulk_structure(element: str, a: float, cubic: bool):
    from ase.build import bulk

    atoms = bulk(name=element, a=a, cubic=cubic)
    # Do I even need to wrap this in a dict?
    return {"structure": atoms}


def generate_structures(
    structure: dict, strain_lst: list[float]
) -> dict[str, dict[str, Atoms]]:
    from ase import Atoms
    structure_lst = []
    structure_ase = Atoms.fromdict(structure)
    for strain in strain_lst:
        structure_strain = structure_ase.copy()
        structure_strain.set_cell(
            structure_strain.cell * strain ** (1 / 3), scale_atoms=True
        )
        structure_lst.append(structure_strain)

    return_dict = {f"qe_{str(i)}": atoms for i, atoms in enumerate(structure_lst)}
    return {"scaled_atoms": return_dict}


def calculate_qe(working_directory, input_dict, structure):
    # ? Possibly don't even define individual functions, but make everything one function

    from pathlib import Path
    from ase import Atoms
    from typing import Any
    import os
    import sys
    import shutil
    import subprocess
    from ase.io import write
    import numpy as np

    def _write_input(input_dict: dict, structure: Atoms, working_directory: str | Path):
        filename = os.path.join(working_directory, "input.pwi")
        os.makedirs(working_directory, exist_ok=True)

        pseudopotentials = input_dict["pseudopotentials"]
        pseudo_path = Path(
            "/home/geiger_j/aiida_projects/adis/git-repos/compare-workflow-graphs/pseudos"
        )
        shutil.copy(src=pseudo_path / pseudopotentials["Al"], dst=working_directory)

        write(
            filename=filename,
            images=structure,
            Crystal=True,
            kpts=input_dict["kpts"],
            input_data={
                "calculation": input_dict["calculation"],
                "occupations": "smearing",
                "degauss": input_dict["smearing"],
                "pseudo_dir": "./",
            },
            pseudopotentials=pseudopotentials,
            tstress=True,
            tprnfor=True,
        )

    def _collect_output(working_directory="."):

        # FIXME: Installed it in OS Python for now, until I know how to use specific Python venv in PythonJob
        from adis_tools.parsers import parse_pw

        output = parse_pw(os.path.join(working_directory, "pwscf.xml"))

        def atoms_to_json_dict(atoms):
            """
            Convert an ASE Atoms object to a fully JSON-serializable dictionary
            that uses only Python base data types.

            Parameters:
            -----------
            atoms : ase.Atoms
                The Atoms object to convert

            Returns:
            --------
            dict
                A dictionary representation using only Python base types
            """
            # Get the dictionary representation from ASE
            atoms_dict = atoms.todict()

            # Create a new dictionary with JSON-serializable values
            json_dict = {}

            # Convert numpy arrays to lists
            for key, value in atoms_dict.items():
                if isinstance(value, np.ndarray):
                    # Convert numpy boolean values to Python booleans
                    if value.dtype == np.bool_ or value.dtype == bool:
                        json_dict[key] = value.tolist()
                    # Convert numpy arrays of numbers to Python lists
                    else:
                        json_dict[key] = value.tolist()
                else:
                    json_dict[key] = value

            return json_dict

        structure_ase = output["ase_structure"]
        structure_json = atoms_to_json_dict(atoms=structure_ase)

        return {
            "structure": structure_json,
            "energy": output["energy"],
            "volume": output["ase_structure"].get_volume(),
        }

    _write_input(
        input_dict=input_dict,
        working_directory=working_directory,
        structure=structure,
    )
    subprocess.check_output(
        "mpirun -np 1 pw.x -in input.pwi > output.pwo",
        cwd=working_directory,
        shell=True,
    )

    outputs = _collect_output(working_directory=working_directory)

    # import ipdb; ipdb.set_trace()

    return outputs


def all_scf(structures, input_dict):
    # Possibly, in this solution, the links of the individual SCF calcs are not resolved in the repr
    from atomistic_engine_unification.calculate_funcs_aiida_wg_pythonjob import calculate_qe

    qe_results = {}
    for key, structure in structures.items():
        # print(key, structure)
        # import ipdb; ipdb.set_trace()
        qe_result = calculate_qe(
            working_directory=key, structure=structure, input_dict=input_dict
        )
        qe_results[key] = qe_result

    return qe_results
    # return {'test': 5}


def plot_energy_volume_curve(qe_results):

    import matplotlib.pyplot as plt

    # import ipdb; ipdb.set_trace()
    energy_lst = [entry['energy'] for entry in qe_results.values()]
    volume_lst = [entry['volume'] for entry in qe_results.values()]

    plt.plot(volume_lst, energy_lst)
    plt.xlabel("Volume")
    plt.ylabel("Energy")
    plt.savefig("evcurve.png")
