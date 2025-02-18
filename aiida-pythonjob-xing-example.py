from aiida_workgraph import WorkGraph, task
from ase.build import bulk
from ase import Atoms
from aiida import load_profile

load_profile()


@task(
    outputs=[
        {"name": "scaled_atoms", "identifier": "workgraph.namespace"},
        {"name": "volumes"},
    ]
)
def generate_scaled_atoms(atoms: Atoms, scales: list) -> dict:
    """Scale the structure by the given scales."""
    volumes = {}
    scaled_atoms = {}
    for i in range(len(scales)):
        atoms1 = atoms.copy()
        atoms1.set_cell(atoms.cell * scales[i], scale_atoms=True)
        scaled_atoms[f"s_{i}"] = atoms1
        volumes[f"s_{i}"] = atoms1.get_volume()
    return {"scaled_atoms": scaled_atoms, "volumes": volumes}


@task()
def emt(atoms):
    from ase.calculators.emt import EMT

    atoms.calc = EMT()
    energy = atoms.get_potential_energy()
    return {"energy": energy}


# Output result from context to the output socket
@task.graph_builder(outputs=[{"name": "results", "from": "context.results"}])
def calculate_enegies(scaled_atoms):
    """Run the scf calculation for each structure."""
    from aiida_workgraph import WorkGraph

    wg = WorkGraph()
    for key, atoms in scaled_atoms.items():
        emt1 = wg.add_task("PythonJob", function=emt, name=f"emt1_{key}", atoms=atoms)
        emt1.set({"computer": "localhost"})
        # save the output parameters to the context
        emt1.set_context({"result": f"results.{key}"})
    return wg


@task()
def fit_eos(volumes: dict, emt_results: dict) -> dict:
    """Fit the EOS of the data."""
    from ase.eos import EquationOfState
    from ase.units import kJ

    volumes_list = []
    energies = []
    for key, data in emt_results.items():
        energy = data["energy"]
        energies.append(energy)
        volumes_list.append(volumes[key])
    #
    eos = EquationOfState(volumes_list, energies)
    v0, e0, B = eos.fit()
    # convert B to GPa
    B = B / kJ * 1.0e24
    eos = {"energy unit": "eV", "v0": v0, "e0": e0, "B": B}
    return eos


atoms = bulk("Au", cubic=True)

wg = WorkGraph("pythonjob_eos_emt")
scale_atoms_task = wg.add_task(
    "PythonJob",
    function=generate_scaled_atoms,
    name="scale_atoms",
    atoms=atoms,
)
# -------- calculate_enegies -----------
calculate_enegies_task = wg.add_task(
    calculate_enegies,
    name="calculate_enegies",
    scaled_atoms=scale_atoms_task.outputs["scaled_atoms"],
)
# -------- fit_eos -----------
wg.add_task(
    "PythonJob",
    function=fit_eos,
    name="fit_eos",
    volumes=scale_atoms_task.outputs["volumes"],
    emt_results=calculate_enegies_task.outputs["results"],
)
# wg.to_html()


wg.submit(
    inputs={
        "scale_atoms": {
            "atoms": atoms,
            "scales": [0.95, 1.0, 1.05],
            "computer": "localhost",
        },
        "fit_eos": {"computer": "localhost"},
    },
    wait=True,
)

print("The fitted EOS parameters are:")
wg.tasks["fit_eos"].outputs["result"].value.value
