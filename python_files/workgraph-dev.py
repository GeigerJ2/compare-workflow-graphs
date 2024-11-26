from aiida_workgraph import WorkGraph

from aiida_workgraph import task


@task
def add(x, y):
    return x + y


@task.calcfunction()
def multiply(x, y):
    return x * y

from aiida import orm, load_profile

load_profile()

wg = WorkGraph("add_multiply_workflow")
wg.add_task(add, name="add1")
wg.add_task(multiply, name="multiply1", x=wg.tasks["add1"].outputs["result"])
# export the workgraph to html file so that it can be visualized in a browser
wg.tasks['add1'].to_html()
wg.to_html()
# visualize the workgraph in jupyter-notebook
# wg
#