from typing import TYPE_CHECKING
import ipdb
from ase import Atoms  # This is only in type hints, so this import should not be necessary

if TYPE_CHECKING:
    from pathlib import Path
    from ase import Atoms


__all__ = ("get_bulk_structure", "calculate_qe")


def get_bulk_structure(element: str, a: float, cubic: bool):
    from ase.build import bulk

    atoms = bulk(name=element, a=a, cubic=cubic)
    # Do I even need to wrap this in a dict?
    return {"structure": atoms}


def generate_structures(
    structure: Atoms, strain_lst: list[float]
) -> dict[str, dict[str, Atoms]]:
    structure_lst = []
    for strain in strain_lst:
        structure_strain = structure.copy()
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

    def _collect_output(working_directory="."):  # , return_atoms=True):
        # FIXME: Installed it in OS Python for now, until I know how to use specific Python venv in PythonJob
        from adis_tools.parsers import parse_pw

        output = parse_pw(os.path.join(working_directory, "pwscf.xml"))
        # if return_atoms:
        return {
            "structure": output["ase_structure"],
            "energy": output["energy"],
            "volume": output["ase_structure"].get_volume(),
        }
        # else:
        #     return {
        #         "energy": output["energy"],
        #         "volume": output["ase_structure"].get_volume(),
        #     }

    # import ipdb; ipdb.set_trace()

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

    import ipdb; ipdb.set_trace()

    return outputs


# ! `to_dict` returns np arrays, and importantly a numpy boolean array for PBC that cannot be serialized:
# ```
# In [3]: atoms
# Out[3]: Atoms(symbols='Al4', pbc=True, cell=[4.05, 4.05, 4.05])

# In [4]: atoms.todict()
# Out[4]:
# {'numbers': array([13, 13, 13, 13]),
# 'positions': array([[0.   , 0.   , 0.   ],
#         [0.   , 2.025, 2.025],
#         [2.025, 0.   , 2.025],
#         [2.025, 2.025, 0.   ]]),
# 'cell': array([[4.05, 0.  , 0.  ],
#         [0.  , 4.05, 0.  ],
#         [0.  , 0.  , 4.05]]),
# 'pbc': array([ True,  True,  True])}

# In [5]: atoms.todict()['pbc'][0]
# Out[5]: True

# In [6]: type(atoms.todict()['pbc'][0])
# Out[6]: numpy.bool_
# ```

# @task.pythonjob()
# ! Possibly pass structure via input_dict
# ! Must be defined inside function that is converted to PythonJob, or discoverable by import from installed module
# def write_input(
#     input_dict: dict[str, Any], structure: Atoms, working_directory: str | Path
# ):
#     filename = os.path.join(working_directory, "input.pwi")
#     os.makedirs(working_directory, exist_ok=True)

#     pseudopotentials = input_dict["pseudopotentials"]
#     shutil.copy(src=pseudo_path / pseudopotentials["Al"], dst=working_directory)

#     write(
#         filename=filename,
#         images=structure,
#         Crystal=True,
#         kpts=input_dict["kpts"],
#         input_data={
#             "calculation": input_dict["calculation"],
#             "occupations": "smearing",
#             "degauss": input_dict["smearing"],
#             "pseudo_dir": "./",
#         },
#         pseudopotentials=pseudopotentials,
#         tstress=True,
#         tprnfor=True,
#     )

# def collect_output(working_directory="."):
#     output = parse_pw(os.path.join(working_directory, "pwscf.xml"))
#     return {
#         "structure": output["ase_structure"].todict(),  # ? `todict`
#         "energy": output["energy"],
#         "volume": output["ase_structure"].get_volume(),
#     }
