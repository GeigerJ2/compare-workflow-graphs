#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pyiron_workflow import Workflow, function_node
from universal_functions import get_kwargs, group_edges_dict, get_source_handles

# In[2]:


from inspect import isfunction


# In[3]:


def create_input_nodes(nodes_dict, edges_lst):
    node_conversion_dict = {
        ed['source']: ed['targetHandle']
        for ed in edges_lst if ed["sourceHandle"] is None
    }
    nodes_to_create_dict = {
        v: nodes_dict[k]
        for k, v in node_conversion_dict.items()
    }
    return nodes_to_create_dict, node_conversion_dict


# In[4]:


def set_input_nodes(workflow, nodes_to_create_dict):
    for k, v in nodes_to_create_dict.items():
        workflow.__setattr__(k, v)
    return workflow


# In[5]:


def get_source_handles(edges_lst):
    source_handle_dict = {}
    for ed in edges_lst:
        if ed['source'] not in source_handle_dict.keys():
            source_handle_dict[ed['source']] = [ed['sourceHandle']]
        else:
            source_handle_dict[ed['source']].append(ed['sourceHandle'])
    return source_handle_dict


# In[6]:


def get_function_nodes(nodes_dict):
    function_dict = {}
    for k, v in nodes_dict.items():
        if isfunction(v):
            if k in source_handle_dict.keys():
                function_dict[k] = {"node_function": v, "output_labels": source_handle_dict[k]}
            else:
                function_dict[k] = {"node_function": v}
    return function_dict


# In[7]:


def get_kwargs(lst):
    return {t['targetHandle']: {'source': t['source'], 'sourceHandle': t['sourceHandle']} for t in lst}




# In[9]:


def build_workflow(workflow, function_dict, total_dict, node_conversion_dict):
    for k, v in function_dict.items():
        kwargs_link_dict = total_dict[k]
        kwargs_dict = {}
        for kw, vw in kwargs_link_dict.items():
            if vw["source"] in node_conversion_dict.keys():
                kwargs_dict[kw] = workflow.__getattribute__(node_conversion_dict[vw["source"]])
            else:
                kwargs_dict[kw] = workflow.__getattr__("tmp_" + str(vw["source"])).__getattribute__("outputs").__getattr__(vw['sourceHandle'])
        v.update(kwargs_dict)
        workflow.__setattr__("tmp_" + str(k), function_node(**v))
    return workflow, "tmp_" + str(k)


# In[10]:


def add_x_and_y(x, y):
    z = x + y
    return x, y, z


# In[11]:


def add_x_and_y_and_z(x, y, z):
    w = x + y + z
    return w


# In[12]:


edges_lst = [
    {'target': 1, 'targetHandle': 'x', 'source': 0, 'sourceHandle': 'x'},
    {'target': 1, 'targetHandle': 'y', 'source': 0, 'sourceHandle': 'y'},
    {'target': 1, 'targetHandle': 'z', 'source': 0, 'sourceHandle': 'z'},
    {'target': 0, 'targetHandle': 'x', 'source': 2, 'sourceHandle': None},
    {'target': 0, 'targetHandle': 'y', 'source': 3, 'sourceHandle': None},
]


# In[13]:


nodes_dict = {
    0: add_x_and_y,
    1: add_x_and_y_and_z,
    2: 1,
    3: 2,
}


# In[14]:


wf = Workflow("my_workflow")


# In[15]:


nodes_to_create_dict, node_conversion_dict = create_input_nodes(nodes_dict=nodes_dict, edges_lst=edges_lst)
wf = set_input_nodes(workflow=wf, nodes_to_create_dict=nodes_to_create_dict)


# In[16]:


source_handle_dict = get_source_handles(edges_lst=edges_lst)
function_dict = get_function_nodes(nodes_dict=nodes_dict)
total_dict = group_edges_dict(edges_lst=edges_lst)


# In[17]:


wf, label = build_workflow(
    workflow=wf,
    function_dict=function_dict,
    total_dict=total_dict,
    node_conversion_dict=node_conversion_dict,
)


# In[18]:


wf.__getattr__(label).pull()


# In[ ]:




