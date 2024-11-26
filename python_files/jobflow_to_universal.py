#!/usr/bin/env python
# coding: utf-8

# In[1]:


from importlib import import_module


# In[2]:


from jobflow import job, Flow
from jobflow.managers.local import run_locally


# In[3]:


@job
def add_xy(x, y):
    return x + y


# In[4]:


@job
def add_xyz(x, y, z):
    return x + y + z


# In[5]:


x = 1
y = 2
z = add_xy(x=x, y=y)
w = add_xyz(x=x, y=y, z=z.output)


# In[6]:


flow = Flow([z, w])


# In[7]:


def get_function_dict(flow):
    return {
        job.uuid: job.function
        for job in flow.jobs
    }


# In[8]:


def get_nodes_dict(function_dict):
    nodes_dict, nodes_mapping_dict = {}, {}
    for i, [k, v] in enumerate(function_dict.items()):
        nodes_dict[i] = v
        nodes_mapping_dict[k] = i
    
    return nodes_dict, nodes_mapping_dict


# In[9]:


def get_edges_and_extend_nodes(flow_dict, nodes_mapping_dict, nodes_dict):
    edges_lst = []
    for job in flow_dict['jobs']:
        for k, v in job['function_kwargs'].items():
            if isinstance(v, dict) and '@module' in v and '@class' in v and '@version' in v:
                edges_lst.append({'target': nodes_mapping_dict[job["uuid"]], 'targetHandle': k, "source": nodes_mapping_dict[v['uuid']], 'sourceHandle': None})
            else:
                if v not in nodes_dict.values():
                    node_index = len(nodes_dict)
                    nodes_dict[node_index] = v
                else:
                    node_index = {tv: tk for tk, tv in nodes_dict.items()}[v]
                edges_lst.append({'target': nodes_mapping_dict[job["uuid"]], 'targetHandle': k, "source": node_index, 'sourceHandle': None})
    return edges_lst, nodes_dict


# In[10]:


flow_dict = flow.as_dict()
function_dict = get_function_dict(flow=flow)
nodes_dict, nodes_mapping_dict = get_nodes_dict(function_dict=function_dict)
edges_lst, nodes_dict = get_edges_and_extend_nodes(flow_dict=flow_dict, nodes_mapping_dict=nodes_mapping_dict, nodes_dict=nodes_dict)
edges_lst, nodes_dict

