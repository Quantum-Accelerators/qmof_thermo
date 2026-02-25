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
)
df = df[df["frac_composition"].apply(lambda d: mol in d)]
df = df.assign(molfrac=df["frac_composition"].apply(lambda d: d[mol]))
df = df[["qmof_id", "pore_diameters", "molfrac", "synthesizable", "chemsys"]]
df["pore_diameters"] = df["pore_diameters"].str[0]

counts = df["synthesizable"].value_counts()
print(counts)

# Separate synthesized and non-synthesized
df_synth = df[df["synthesizable"] == True]
df_nonsynth = df[df["synthesizable"] == False]

fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(df_synth["molfrac"], df_synth["pore_diameters"],
           c='#6495ED', alpha=1, s=20, label=f'Synthesized ({len(df_synth)})')

# Plot non-synthesized first (background)
ax.scatter(df_nonsynth["molfrac"], df_nonsynth["pore_diameters"],
           facecolors='none', edgecolors='#FF6347', alpha=1, s=20, linewidths=0.4, label=f'Not Synthesized ({len(df_nonsynth)})')

ax.text(-0.145, 1.0, "B", transform=ax.transAxes, fontsize=20,
         fontweight="bold", va="top", ha="left")


ax.tick_params(which='major', direction='in', length=26, width=1.8)
ax.tick_params(which='minor', direction='in', length=14, width=1.6)
ax.tick_params(labelsize=20)
ax.minorticks_on()

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

plt.xlabel("Carbon Fraction in MOF", fontsize=22)
plt.ylabel("Largest Cavity Diameter (Å)", fontsize=22)
plt.xlim([0, 0.7])
plt.ylim([0, 48])
plt.legend(fontsize=16, frameon=True, loc='upper right')

plt.tight_layout()
plt.savefig("FigureS25b.png", bbox_inches='tight', dpi=500)
plt.show()
