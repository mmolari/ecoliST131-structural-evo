# %%
import pandas as pd
import subprocess
import pathlib
from Bio import Phylo
import matplotlib.pyplot as plt


dset = "ST131_ABC"


fig_fld = pathlib.Path("figs/f06")
fig_fld.mkdir(exist_ok=True, parents=True)


# load joint coordinates dictionary
dset = "ST131_ABC"
fld = pathlib.Path(f"../../results/{dset}")

# %%
junction = "ANHYQUDQAA_r__IMAMHFLCWS_f"
# junction = "ATFHQYFPNW_f__GUDQDOMJFJ_f"
# junction = "JVNRLCFAVD_f__PLTCZQCVRD_r"


pan_file = fld / f"backbone_joints/asm20-100-5/joints_pangraph/{junction}.json"

exp_fld = fig_fld / junction

# pangraph export
command = f"pangraph export {pan_file} -o {exp_fld} -ell 0 -nd -ll 0"
subprocess.run(command, shell=True, check=True)

# %%

# load tree
tree_file = fld / "pangraph/asm20-100-5-filtered-coretree.nwk"
tree = Phylo.read(tree_file, "newick")
tree.ladderize()

# draw tree
fig, ax = plt.subplots(figsize=(3, 10))
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: "")
ax.set_axis_off()
ax.plot([0, 1e-5], [100, 100], color="k", linewidth=1)
ax.text(0.5e-5, 99, "1e-5", ha="center", va="bottom", fontsize=10)

fig.savefig(fig_fld / "tree.svg", bbox_inches="tight", pad_inches=0)
plt.show()

# %%
