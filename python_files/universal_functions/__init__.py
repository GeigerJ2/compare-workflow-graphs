from importlib import import_module
from inspect import isfunction
from rich.pretty import pprint
from jobflow import job
from aiida_workgraph import task

__all__ = (
    "get_kwargs",
    "group_edges",
    "get_input_dict",
    "get_workflow",
    "get_source_handles",
    "remove_source_handles",
    "universal_reduce",
)


def get_kwargs(lst):
    return {
        t["targetHandle"]: {"source": t["source"], "sourceHandle": t["sourceHandle"]}
        for t in lst
    }


def group_edges_dict(edges_lst):
    # edges_sorted_lst = sorted(edges_lst, key=lambda x: x['target'], reverse=True)
    total_dict = {}
    tmp_lst = []
    target_id = edges_lst[0]["target"]
    for ed in edges_lst:
        if target_id == ed["target"]:
            tmp_lst.append(ed)
        else:
            total_dict[target_id] = get_kwargs(lst=tmp_lst)
            target_id = ed["target"]
            tmp_lst = [ed]
    total_dict[target_id] = get_kwargs(lst=tmp_lst)
    return total_dict


def group_edges_list(edges_lst):
    edges_sorted_lst = sorted(edges_lst, key=lambda x: x["target"], reverse=True)
    total_lst, tmp_lst = [], []
    target_id = edges_sorted_lst[0]["target"]
    for ed in edges_sorted_lst:
        if target_id == ed["target"]:
            tmp_lst.append(ed)
        else:
            total_lst.append((target_id, get_kwargs(lst=tmp_lst)))
            target_id = ed["target"]
            tmp_lst = [ed]
    total_lst.append((target_id, get_kwargs(lst=tmp_lst)))
    return total_lst


def get_input_dict(nodes_dict):
    return {k: v for k, v in nodes_dict.items() if not isfunction(v)}



def get_source_handles(edges_lst):
    source_handle_dict = {}
    for ed in edges_lst:
        if ed["source"] not in source_handle_dict.keys():
            source_handle_dict[ed["source"]] = [ed["sourceHandle"]]
        else:
            source_handle_dict[ed["source"]].append(ed["sourceHandle"])
    return source_handle_dict


def remove_source_handles(nodes_dict, edges_lst):
    nodes_updated_dict = {}
    translate_dict = {}
    i = 0
    for k, v in nodes_dict.items():
        if not isfunction(v):
            nodes_updated_dict[i] = v
            translate_dict[k] = i
            i += 1

    edges_updated_lst = []
    i = len(nodes_updated_dict)
    source_handle_dict = get_source_handles(edges_lst=edges_lst)
    for k, v in nodes_dict.items():
        if k not in translate_dict.keys():
            edges_to_be_updated_dict = {}
            for ed in edges_lst:
                if ed["target"] == k and ed["sourceHandle"] is None:
                    edges_updated_lst.append(
                        {
                            "target": i,
                            "targetHandle": ed["targetHandle"],
                            "source": translate_dict[ed["source"]],
                            "sourceHandle": None,
                        }
                    )
                elif ed["target"] == k and ed["sourceHandle"] is not None:
                    nodes_updated_dict[i] = ed["sourceHandle"]
                    source_handle_index = i
                    i += 1
                    new_item = source_handle_dict[ed["source"]]
                    if new_item not in nodes_updated_dict.values():
                        nodes_updated_dict[i] = new_item
                        output_label_index = i
                        i += 1
                    else:
                        for kn, vn in nodes_updated_dict.items():
                            if vn == new_item:
                                output_label_index = kn
                    nodes_updated_dict[i] = get_item_from_tuple
                    edges_updated_lst.append(
                        {
                            "target": i,
                            "targetHandle": "input_obj",
                            "source": translate_dict[ed["source"]],
                            "sourceHandle": None,
                        }
                    )
                    edges_updated_lst.append(
                        {
                            "target": i,
                            "targetHandle": "index",
                            "source": source_handle_index,
                            "sourceHandle": None,
                        }
                    )
                    edges_updated_lst.append(
                        {
                            "target": i,
                            "targetHandle": "index_lst",
                            "source": output_label_index,
                            "sourceHandle": None,
                        }
                    )
                    edges_to_be_updated_dict[i] = ed
                    i += 1
            nodes_updated_dict[i] = v
            translate_dict[k] = i
            for k, v in edges_to_be_updated_dict.items():
                edges_updated_lst.append(
                    {
                        "target": i,
                        "targetHandle": v["targetHandle"],
                        "source": k,
                        "sourceHandle": None,
                    }
                )
            i += 1

    return nodes_updated_dict, edges_updated_lst


def universal_reduce(nodes_dict, edges_lst):
    edges_with_source_handles_lst = [
        ed for ed in edges_lst if ed["sourceHandle"] is not None
    ]
    if len(edges_with_source_handles_lst) > 0:
        return remove_source_handles(nodes_dict=nodes_dict, edges_lst=edges_lst)
    else:
        return nodes_dict, edges_lst


# @task
def get_item_from_tuple(input_obj, index, index_lst):
    if isinstance(input_obj, dict):
        return input_obj[index]
    else:  # input_obj is a tuple
        return list(input_obj)[index_lst.index(index)]



def get_workflow(nodes_dict, input_dict, total_dict, backend: str = 'aiida'):
    memory_dict = {}

    pprint(nodes_dict)
    pprint(input_dict)
    pprint(total_dict)

    for k, v in nodes_dict.items():
        print(f"k: {k}")
        pprint(v)
        print("*" * 40)
        if isfunction(v):
            if backend == 'jobflow':
                fn = job(v)
                output_identifier = 'output'
            elif backend == 'aiida':
                fn = task(v)
                output_identifier = 'result'
            else:
                raise SystemExit

            kwargs = {}
            # print(f"total_dict[k]: {total_dict[k]}")
            for kw, vw in total_dict[k].items():
                print(f"kw: {kw}")
                pprint(vw)
                print("#" * 40)
                if vw["source"] in input_dict:
                    kwargs[kw] = input_dict[vw["source"]]
                else:
                    kwargs[kw] = getattr(memory_dict[vw["source"]], output_identifier)
            # kwargs = {
            #     kw: input_dict[vw["source"]]
            #     if vw["source"] in input_dict
            #     else getattr(memory_dict[vw["source"]], output_identifier)
            #     for kw, vw in total_dict[k].items()
            # }
            memory_dict[k] = fn(**kwargs)
    return list(memory_dict.values())
