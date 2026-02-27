from __future__ import annotations

import gzip
import json

import matplotlib.pyplot as plt
import pandas as pd

# 1) Load your results JSON
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# 2) Build a tidy DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
    .reset_index()
    .rename(columns={"index": "qmof_id"})
)

df = df[["ehull", "pore_diameters", "synthesizable", "chemsys"]]
# df = df[[ "ehull", "pore_diameters", "synthesizable"]]
df["pore_diameters"] = df["pore_diameters"].str[0]

counts = df["synthesizable"].value_counts()
print(counts)

df = df[df["synthesizable"]]

# df = df[df["chemsys"] == "C-H-N-Zn"]

print(df)


counts = df["synthesizable"].value_counts()
print(counts)


color_map = {True: "#6495ED", False: "#FF6347"}
colors = df["synthesizable"].map(color_map)

# 5) Create scatter plot
fig, ax = plt.subplots(figsize=(3.4, 2.485))
plt.scatter(
    df["pore_diameters"], df["ehull"], c=colors, s=4, alpha=0.4, edgecolors="none"
)
plt.xlabel("Largest Cavity Diameter (Å)", fontsize=8)
plt.ylabel(r"$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=8)


ax.tick_params(which="major", direction="in", length=10, width=1.25)
ax.tick_params(which="minor", direction="in", length=5, width=1.25)

for spine in ax.spines.values():
    spine.set_linewidth(1.25)

ax.tick_params(labelsize=8)

ax.minorticks_on()

# plt.legend(handles=[
#    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#6495ED', markersize=12, label='Synthesized'),
#    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6347', markersize=12, label='Not Synthesized')
# ],
# fontsize = 22,
# loc = 'best',
# frameon = False
# )
plt.xlim([0, 38])
plt.ylim([0, 0.82])

# ax.text(
#    -0.145,    # x position, just left of the left spine
#    1.0,     # y position, just above the top spine
#    "B",      # the label
#    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
#    fontsize=11,              # size of the letter
#    fontweight="bold",        # make it bold
#    va="top",                 # vertical alignment
#    ha="left"                 # horizontal alignment
# )

# plt.xlim([0, 40])
# plt.ylim([0, 0.82])

plt.tight_layout()
plt.savefig("FigureS15", dpi=1500, bbox_inches="tight")
plt.show()
# s#ns.set_theme(style="ticks")
