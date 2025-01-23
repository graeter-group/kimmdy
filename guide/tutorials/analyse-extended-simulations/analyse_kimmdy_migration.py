# %%
import MDAnalysis as mda
from pymol import cmd
from pymol import cgo
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
import json


# %%
def get_inbetween(p1, p2, val):
    return p1 + val * (p2 - p1)


# %%
interactions_path = Path("radical_migration_dopa_merged.json")
geometry_path = Path("dopa.gro")
do_arrow = True
analysis_type = "count"
manual_max_count = 40

# %%
u = mda.Universe(geometry_path.as_posix(), guess_bonds=True)

with open(interactions_path, "r") as json_file:
    interactions = json.load(json_file)

# %%
rgb = None
if analysis_type == "quantitative":
    rgb = np.tile(np.asarray([0.62, 0.59, 0.98]), (len(interactions.keys()), 1))
elif analysis_type == "max_rate":
    max_rates = np.array([entry["max_rate"] for entry in interactions.values()])
    log_norm_max_rates = np.log(max_rates)

    # Normalize the log values to the range [0, 1]
    min_log = np.min(log_norm_max_rates)
    max_log = np.max(log_norm_max_rates)
    norm_log_max_rates = (log_norm_max_rates - min_log) / (max_log - min_log)

    # Use the purple sequential colormap
    cmap = plt.get_cmap("Purples")

    # Get RGB values for the normalized log max_rates
    rgb = [np.asarray(cmap(value)[:3]) for value in norm_log_max_rates]
elif analysis_type == "count":
    print("Selected analysis by count of reaction occurence.")

    counts = np.array([entry["count"] for entry in interactions.values()])

    cmap = mpl.colormaps["autumn"].reversed()
    norm = mpl.colors.Normalize(vmin=counts.min(), vmax=counts.max())
    rgb = [np.asarray(cmap(norm(value))[:3]) for value in counts]


# %%
radius = 10
counter = 0
for k, v in interactions.items():
    print(f"Building object for atoms {k}")
    atom_1_nr, atom_2_nr = k.split(sep="_")

    a1 = u.select_atoms(f"id {atom_1_nr}")[0].position
    a2 = u.select_atoms(f"id {atom_2_nr}")[0].position

    p_cap = get_inbetween(a1, a2, 0.75)
    radius = 0.12
    p1 = get_inbetween(a1, a2, 0.1)
    p2 = get_inbetween(a1, a2, 0.9)

    if do_arrow and rgb:
        p1_far = get_inbetween(a1, a2, 0.18)

        cylinder = (
            [cgo.CYLINDER]
            + p1_far.tolist()
            + p_cap.tolist()
            + [radius]
            + rgb[counter].tolist()
            + rgb[counter].tolist()
        )
        head = (
            [cgo.CONE]
            + p_cap.tolist()
            + p2.tolist()
            + [radius * 2, 0.0]
            + rgb[counter].tolist()
            + rgb[counter].tolist()
            + [1.0, 0.0]
        )
        # cmd.load_cgo(cylinder,f"cylinder{i}")
        cmd.load_cgo(cylinder + head, f"cylinder{k}")
        counter += 1
print(f"Loading geometry file")
cmd.load(geometry_path.as_posix())
print("Changing default settings")
cmd.set("virtual_trackball", 0)
cmd.bg_color("white")
cmd.set("orthoscopic", 1)
cmd.set("stick_radius", 0.1)
cmd.show("licorice", "all")
cmd.hide("everything", "resn SOL")
# %%
