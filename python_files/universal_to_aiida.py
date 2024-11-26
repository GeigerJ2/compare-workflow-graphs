# In[ ]:


from aiida import orm, load_profile

load_profile()

from aiida_workgraph import task, WorkGraph
from jobflow import job, Flow
from jobflow.managers.local import run_locally

from rich.pretty import pprint

from universal_functions import group_edges_dict, get_input_dict, get_workflow, universal_reduce

# In[ ]:

# Actual tasks


# @task.calcfunction()
def add_xy(x, y):
    return x + y, x, y


# @task.calcfunction()
def add_xyz(x, y, z):
    return x + y + z




# In[ ]:

# Actually running the workflow

nodes_dict = {
    0: add_xy,
    1: add_xyz,
    2: 1,
    3: 2,
}


edges_lst = [
    {"target": 1, "targetHandle": "x", "source": 0, "sourceHandle": "x"},
    {"target": 1, "targetHandle": "y", "source": 0, "sourceHandle": "y"},
    {"target": 1, "targetHandle": "z", "source": 0, "sourceHandle": "z"},
    {"target": 0, "targetHandle": "x", "source": 2, "sourceHandle": None},
    {"target": 0, "targetHandle": "y", "source": 3, "sourceHandle": None},
]


nodes_updated_dict, edges_updated_lst = universal_reduce(
    nodes_dict=nodes_dict, edges_lst=edges_lst
)
nodes_updated_dict, edges_updated_lst


total_dict = group_edges_dict(edges_lst=edges_updated_lst)
input_dict = get_input_dict(nodes_dict=nodes_updated_dict)

pprint(total_dict)
# raise SystemExit

task_lst = get_workflow(
    nodes_dict=nodes_updated_dict,
    input_dict=input_dict,
    total_dict=total_dict,
    backend='aiida'
)

raise SystemExit


flow = Flow(task_lst)


# In[ ]:


result = run_locally(flow)
result


# In[ ]:
