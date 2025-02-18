from pyiron_base import job

__all__ = (
    "get_delayed_object_dict",
    "get_dict",
    "get_kwargs",
    "get_list",
    "get_source",
    "get_source_handles",
    "group_edges",
    "resort_total_lst",
)


def get_delayed_object_dict(total_lst, nodes_dict, source_handle_dict, pyiron_project):
    delayed_object_dict = {}
    for item in total_lst:
        key, input_dict = item
        kwargs = {
            k: get_source(
                nodes_dict=nodes_dict,
                delayed_object_dict=delayed_object_dict,
                source=v["source"],
                sourceHandle=v["sourceHandle"],
            )
            for k, v in input_dict.items()
        }
        # print(nodes_dict[key], source_handle_dict.get(key, []))
        # print(kwargs)
        delayed_object_dict[key] = job(
            funct=nodes_dict[key],
            output_key_lst=source_handle_dict.get(key, []),
        )(**kwargs, pyiron_project=pyiron_project)
    return delayed_object_dict


def get_kwargs(lst):
    return {
        t["targetHandle"]: {"source": t["source"], "sourceHandle": t["sourceHandle"]}
        for t in lst
    }


def resort_total_lst(total_lst, nodes_dict):
    nodes_with_dep_lst = list(sorted([v[0] for v in total_lst]))
    nodes_without_dep_lst = [
        k for k in nodes_dict.keys() if k not in nodes_with_dep_lst
    ]
    ordered_lst, total_new_lst = [], []
    while len(total_new_lst) < len(total_lst):
        for ind, connect in total_lst:
            if ind not in ordered_lst:
                source_lst = [sd["source"] for sd in connect.values()]
                if all(
                    [s in ordered_lst or s in nodes_without_dep_lst for s in source_lst]
                ):
                    ordered_lst.append(ind)
                    total_new_lst.append([ind, connect])
    return total_new_lst


def group_edges(edges_lst):
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


def get_source_handles(edges_lst):
    source_handle_dict = {}
    for ed in edges_lst:
        if ed["source"] not in source_handle_dict.keys():
            source_handle_dict[ed["source"]] = [ed["sourceHandle"]]
        else:
            source_handle_dict[ed["source"]].append(ed["sourceHandle"])
    return {
        k: list(range(len(v))) if len(v) > 1 and all([el is None for el in v]) else v
        for k, v in source_handle_dict.items()
    }


def get_source(nodes_dict, delayed_object_dict, source, sourceHandle):
    if source in delayed_object_dict.keys():
        return (
            delayed_object_dict[source].__getattr__("output").__getattr__(sourceHandle)
        )
    else:
        return nodes_dict[source]


def get_dict(**kwargs):
    return {k: v for k, v in kwargs["kwargs"].items()}


def get_list(**kwargs):
    return list(kwargs["kwargs"].values())
