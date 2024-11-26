#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyiron_base
pyiron_base.__file__


# In[2]:


from pyiron_base import Project, job
from pyiron_base.project.delayed import draw, DelayedObject


# In[3]:


def remove_server_obj(nodes_dict, edges_lst):
    server_lst = [k for k in nodes_dict.keys() if k.startswith("server_obj_")]
    for s in server_lst:
        del nodes_dict[s]
        edges_lst = [ep for ep in edges_lst if s not in ep]
    return nodes_dict, edges_lst


# In[4]:


def get_nodes(connection_dict, delayed_object_updated_dict):
    return {connection_dict[k]: v._python_function if isinstance(v, DelayedObject) else v for k, v in delayed_object_updated_dict.items()}


# In[5]:


def get_unique_objects(nodes_dict, edges_lst):
    delayed_object_dict = {k:v for k, v in nodes_dict.items() if isinstance(v, DelayedObject)}
    unique_lst = []
    delayed_object_updated_dict, match_dict = {}, {}
    for dobj in delayed_object_dict.keys():
        match = False
        for obj in unique_lst:
            if delayed_object_dict[dobj]._input == delayed_object_dict[obj]._input:
                delayed_object_updated_dict[obj] = delayed_object_dict[obj]
                match_dict[dobj] = obj
                match = True
                break
        if not match:
            unique_lst.append(dobj)
            delayed_object_updated_dict[dobj] = delayed_object_dict[dobj]
    delayed_object_updated_dict.update({k:v for k, v in nodes_dict.items() if not isinstance(v, DelayedObject)})
    return delayed_object_updated_dict, match_dict


# In[6]:


def get_connection_dict(delayed_object_updated_dict, match_dict):
    new_obj_dict = {}
    connection_dict = {}
    lookup_dict = {}
    for i, [k, v] in enumerate(delayed_object_updated_dict.items()):
        new_obj_dict[i] = v
        connection_dict[k] = i
        lookup_dict[i] = k

    for k, v in match_dict.items():
        if v in connection_dict.keys():
            connection_dict[k] = connection_dict[v]

    return connection_dict, lookup_dict


# In[7]:



# # Define Workflow

# In[8]:


pr = Project("test")


# In[9]:


pr.remove_jobs(recursive=True, silently=True)


# In[10]:


@job(output_key_lst=["x", "y", "z"])
def add_x_and_y(x, y):
    z = x + y
    return {"x": x, "y": y, "z": z}


# In[11]:


@job
def add_x_and_y_and_z(x, y, z):
    w = x + y + z
    return w


# In[12]:


pr.remove_jobs(recursive=True, silently=True)


# In[13]:


obj = add_x_and_y(x=1, y=2, pyiron_project=pr)


# In[14]:


w = add_x_and_y_and_z(x=obj.output.x, y=obj.output.y, z=obj.output.z, pyiron_project=pr)
w.pull()


# In[15]:


pr.job_table()


# In[16]:


w = add_x_and_y_and_z(x=obj.output.x, y=obj.output.y, z=obj.output.z, pyiron_project=pr)
w.pull()


# # Convert to universal format

# In[17]:


nodes_dict, edges_lst = w.get_graph()


# In[18]:


nodes_dict, edges_lst = remove_server_obj(nodes_dict=nodes_dict, edges_lst=edges_lst)


# In[19]:


delayed_object_updated_dict, match_dict = get_unique_objects(nodes_dict=nodes_dict, edges_lst=edges_lst)


# In[20]:


connection_dict, lookup_dict = get_connection_dict(delayed_object_updated_dict=delayed_object_updated_dict, match_dict=match_dict)


# In[21]:


get_nodes(connection_dict=connection_dict, delayed_object_updated_dict=delayed_object_updated_dict)


# In[22]:


get_edges_dict(edges_lst=edges_lst, nodes_dict=nodes_dict, connection_dict=connection_dict, lookup_dict=lookup_dict)


# In[ ]:




