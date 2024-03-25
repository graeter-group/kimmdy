from kimmdy.config import Config
from kimmdy.runmanager import RunManager
from pathlib import Path

# reimplement this without checkpoints
# def test_tasks_are_set_up(arranged_tmp_path):
#     """
#     use the initial checkpoint writing
#     to test properties of the runmanager
#     initialization without having to `.run()` it.
#     """
#     config = Config(Path("config1.yml"))
#     runmgr = RunManager(config)
#     runmgr.write_one_checkpoint()

#     items = []
#     while not runmgr.tasks.empty():
#         items.append(runmgr.tasks.get().name)
#     assert items == [
#         "_setup",
#         "_run_md",
#         "_run_md",
#         "_place_reaction_tasks",
#         "_decide_recipe",
#         "_apply_recipe",
#     ]

#     config2 = Config(Path("config2.yml"))
#     runmgr2 = RunManager(config2)
#     runmgr2.write_one_checkpoint()

#     items = []
#     while not runmgr2.tasks.empty():
#         items.append(runmgr2.tasks.get().name)
#     assert items == [
#         "_setup",
#         "_run_md",
#         "_run_md",
#         "_place_reaction_tasks",
#         "_decide_recipe",
#         "_apply_recipe",
#         "_place_reaction_tasks",
#         "_decide_recipe",
#         "_apply_recipe",
#         "_run_md",
#     ]
