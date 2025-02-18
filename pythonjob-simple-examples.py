from aiida import load_profile
from aiida.engine import run_get_node
from aiida_pythonjob import PythonJob, prepare_pythonjob_inputs

load_profile()


def add(x, y):
    return x + y


inputs = prepare_pythonjob_inputs(
    add,
    function_inputs={"x": 1, "y": 2},
    computer="localhost",
)
result, node = run_get_node(PythonJob, inputs=inputs)
print("result: ", result["result"])