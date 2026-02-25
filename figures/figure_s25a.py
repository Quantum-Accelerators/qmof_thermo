import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import gzip

with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

mol = "C"
# build base df and extract co2_frac
df = (
    pd.DataFrame.from_dict(data, orient="index")
      .rename_axis("qmof_id")
#      .reset_index()
)
df = df[df["frac_composition"].apply(lambda d: mol in d)]
df = df.assign(molfrac=df["frac_composition"].apply(lambda d: d[mol]))
df = df[["qmof_id", "pore_diameters", "molfrac", "synthesizable", "chemsys"]]
#df = df[df["chemsys"]  == "C-H-N-O-Zn"]
df["pore_diameters"] = df["pore_diameters"].str[0]


counts = df["synthesizable"].value_counts()
print(counts)

fig, ax = plt.subplots(figsize=(8, 6))

# Create hexbin plot with log-scale colorbar
hb = ax.hexbin(df["molfrac"], df["pore_diameters"],
               gridsize=50, cmap='viridis', mincnt=1,
               extent=([0, 0.7, 0, 48]),
               norm=LogNorm())

# Add colorbar
cb = plt.colorbar(hb, ax=ax, pad=0.02)

cb.ax.tick_params(labelsize=20)


ax.tick_params(
    which='major', direction='in', length=26, width=1.8
)
ax.tick_params(
    which='minor', direction='in', length=14, width=1.6
)
ax.tick_params(labelsize=20)
ax.minorticks_on()

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

plt.xlabel("Carbon Fraction in MOF", fontsize=22)
plt.ylabel("Largest Cavity Diameter (Å)", fontsize=22)

ax.text(-0.145, 1.0, "A", transform=ax.transAxes, fontsize=20,
         fontweight="bold", va="top", ha="left")


#plt.xlim([0, 0.68])
#plt.ylim([-0.06, 0.82])

plt.tight_layout()  # make sure things are laid out nicely
plt.savefig("FigureS25a.png", bbox_inches='tight', dpi=500)
plt.show()

