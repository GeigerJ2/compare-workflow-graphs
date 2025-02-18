# FIXME: These are still all the pyiron functions

# from pyiron_base.project.delayed import DelayedObject

# __all__ = (
#     'get_connection_dict',
#     'get_dict',
#     'get_edges_dict',
#     'get_list',
#     'get_nodes',
#     'get_unique_objects',
#     'remove_server_obj',
# )


# def remove_server_obj(nodes_dict, edges_lst):
#     server_lst = [k for k in nodes_dict.keys() if k.startswith("server_obj_")]
#     for s in server_lst:
#         del nodes_dict[s]
#         edges_lst = [ep for ep in edges_lst if s not in ep]
#     return nodes_dict, edges_lst


# def get_nodes(connection_dict, delayed_object_updated_dict):
#     return {
#         connection_dict[k]: v._python_function if isinstance(v, DelayedObject) else v
#         for k, v in delayed_object_updated_dict.items()
#     }


# def get_unique_objects(
#     nodes_dict, edges_lst
# ):  # I need a pre-filter before this function to remove the virtual nodes
#     delayed_object_dict = {}
#     for k, v in nodes_dict.items():
#         if isinstance(v, DelayedObject):
#             delayed_object_dict[k] = v
#         elif (
#             isinstance(v, list) and any([isinstance(el, DelayedObject) for el in v])
#         ):  # currently this replaces any list - what I need instead is some kind of virtual node - mixed nodes
#             delayed_object_dict[k] = DelayedObject(function=get_list)
#             delayed_object_dict[k]._input = {i: el for el in enumerate(v)}
#             delayed_object_dict[k]._python_function = get_list
#         elif isinstance(v, dict) and any(
#             [isinstance(el, DelayedObject) for el in v.values()]
#         ):
#             delayed_object_dict[k] = DelayedObject(
#                 function=get_dict,
#                 **v,
#             )
#             delayed_object_dict[k]._python_function = get_dict
#             delayed_object_dict[k]._input = v
#     unique_lst = []
#     delayed_object_updated_dict, match_dict = {}, {}
#     for dobj in delayed_object_dict.keys():
#         match = False
#         for obj in unique_lst:
#             # print(delayed_object_dict[dobj]._list_index, delayed_object_dict[dobj]._output_key, delayed_object_dict[obj]._list_index, delayed_object_dict[obj]._output_key)
#             if (
#                 obj.split("_")[0] == dobj.split("_")[0]
#                 and delayed_object_dict[dobj]._input == delayed_object_dict[obj]._input
#             ):
#                 delayed_object_updated_dict[obj] = delayed_object_dict[obj]
#                 match_dict[dobj] = obj
#                 match = True
#                 break
#         if not match:
#             unique_lst.append(dobj)
#             delayed_object_updated_dict[dobj] = delayed_object_dict[dobj]
#     update_dict = {}
#     for k, v in nodes_dict.items():
#         if not (
#             isinstance(v, DelayedObject)
#             or (
#                 isinstance(v, list) and any([isinstance(el, DelayedObject) for el in v])
#             )
#             or (
#                 isinstance(v, dict)
#                 and any([isinstance(el, DelayedObject) for el in v.values()])
#             )
#         ):
#             update_dict[k] = v
#     delayed_object_updated_dict.update(update_dict)
#     return delayed_object_updated_dict, match_dict


# def get_connection_dict(delayed_object_updated_dict, match_dict):
#     new_obj_dict = {}
#     connection_dict = {}
#     lookup_dict = {}
#     for i, [k, v] in enumerate(delayed_object_updated_dict.items()):
#         new_obj_dict[i] = v
#         connection_dict[k] = i
#         lookup_dict[i] = k

#     for k, v in match_dict.items():
#         if v in connection_dict.keys():
#             connection_dict[k] = connection_dict[v]

#     return connection_dict, lookup_dict


# def get_edges_dict(edges_lst, nodes_dict, connection_dict, lookup_dict):
#     edges_dict_lst = []
#     existing_connection_lst = []
#     for ep in edges_lst:
#         input_name, output_name = ep
#         target = connection_dict[input_name]
#         target_handle = "_".join(output_name.split("_")[:-1])
#         connection_name = lookup_dict[target] + "_" + target_handle
#         if connection_name not in existing_connection_lst:
#             output = nodes_dict[output_name]
#             if isinstance(output, DelayedObject):
#                 if output._list_index is not None:
#                     edges_dict_lst.append(
#                         {
#                             "target": target,
#                             "targetHandle": target_handle,
#                             "source": connection_dict[output_name],
#                             "sourceHandle": output._list_index,  # check for list index
#                         }
#                     )
#                 else:
#                     edges_dict_lst.append(
#                         {
#                             "target": target,
#                             "targetHandle": target_handle,
#                             "source": connection_dict[output_name],
#                             "sourceHandle": output._output_key,  # check for list index
#                         }
#                     )
#             else:
#                 edges_dict_lst.append(
#                     {
#                         "target": target,
#                         "targetHandle": target_handle,
#                         "source": connection_dict[output_name],
#                         "sourceHandle": None,
#                     }
#                 )
#             existing_connection_lst.append(connection_name)
#     return edges_dict_lst


# def get_dict(**kwargs):
#     return {k: v for k, v in kwargs["kwargs"].items()}


# def get_list(**kwargs):
#     return list(kwargs["kwargs"].values())
