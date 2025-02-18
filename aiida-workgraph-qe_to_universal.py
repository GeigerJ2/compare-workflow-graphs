import ipdb
from ase.build import bulk
import pathlib

from aiida import load_profile, orm
from aiida_workgraph import WorkGraph, task

from to_universal_funcs_workgraph_calcfunctions import (
    calculate_qe,
    collect_output,
    generate_structures,
    get_bulk_structure,
    write_input,
    run_qe,
    all_scf,
    plot_energy_volume_curve,
)

load_profile()

pseudopotentials = {"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"}
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


# ? Why does
# get_bulk_structure = task.calcfunction(outputs=[{"name": "structure"}])(
#     get_bulk_structure
# )

@task.calcfunction(outputs=[{"name": "structure"}])
def get_bulk_structure(element: orm.Str, a: orm.Float, cubic: orm.Bool):
    atoms = bulk(name=element.value, a=a.value, cubic=cubic.value)
    return {"structure": orm.StructureData(ase=atoms)}

generate_structures = task.calcfunction(outputs=[{"name": "structures"}])(
    generate_structures
)


collect_output = task.calcfunction(
    outputs=[{"name": "structure"}, {"name": "energy"}, {"name": "volume"}]
)(collect_output)


write_input = task.calcfunction()(write_input)


run_qe = task.calcfunction()(run_qe)


# ! Or, maybe, just make the whole process of writing the input and running the QE process as one calcfunction
calculate_qe = task.graph_builder(
    outputs=[
        {"name": "structure", "from": "collect_output.structure"},
        {"name": "energy", "from": "collect_output.energy"},
        {"name": "volume", "from": "collect_output.volume"},
    ]
)(calculate_qe)


all_scf = task.graph_builder(
    outputs=[{"name": "result", "from": "context.result"}],
)(all_scf)


plot_energy_volume_curve = task.calcfunction()(plot_energy_volume_curve)


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

# ipdb.set_trace()
wg_eos.run()

raise SystemExit()

# ? Here, convert to universal representation

work_graph_dict = wg_eos.to_dict()

# ! This seems eerily similar to the `NodeLink` class from node-graph
edges_label_lst = []
for link_dict in work_graph_dict["links"]:
    if link_dict["from_socket"] == "result":
        edges_label_lst.append(
            {
                "target": link_dict["to_node"],
                "targetHandle": link_dict["to_socket"],
                "source": link_dict["from_node"],
                "sourceHandle": None,
            }
        )
    else:
        edges_label_lst.append(
            {
                "target": link_dict["to_node"],
                "targetHandle": link_dict["to_socket"],
                "source": link_dict["from_node"],
                "sourceHandle": link_dict["from_socket"],
            }
        )

print(edges_label_lst)

kwargs_dict, function_dict = {}, {}
for task_name, task_dict in work_graph_dict["tasks"].items():
    # try:
    #     input_variables = [
    #         input_parameter
    #         for input_parameter in task_dict['inputs'].keys()
    #         if not input_parameter.startswith("metadata") and not input_parameter.startswith("_wait")
    #     ]
    #     input_kwargs = {
    #         input_parameter: task_dict['inputs'][input_parameter]['property']["value"].value if isinstance(task_dict['inputs'][input_parameter]['property']["value"], dict) else task_dict['inputs'][input_parameter]['property']["value"]
    #         for input_parameter in input_variables
    #     }
    #     kwargs_dict[task_name] = input_kwargs

    # Extract relevant input variable names (excluding metadata and _wait)
    input_variables = [
        param
        for param in task_dict["inputs"]
        if not param.startswith("metadata") and not param.startswith("_wait")
    ]

    # Prepare input keyword arguments
    input_kwargs = {}

    for param in input_variables:
        print(f"PARAM: {param}")
        try:
            property_value = task_dict["inputs"][param]["property"]["value"]
        except:
            ipdb.set_trace()
            raise

        if isinstance(property_value, dict):
            input_kwargs[param] = (
                property_value.value
            )  # Assuming `.value` is a valid key
        else:
            input_kwargs[param] = property_value

    # function_dict[task_name] = loads(task_dict['executor']['callable']).process_class._func
    print(f"{task_name.upper()}: {task_dict['executor']['callable']}")
    kwargs_dict[task_name] = input_kwargs
    function_dict[task_name] = task_dict["executor"]["callable"]


from rich.pretty import pprint

print("KWARGS_DICT")
pprint(kwargs_dict)
print("FUNCTION_DICT")
print(function_dict)
