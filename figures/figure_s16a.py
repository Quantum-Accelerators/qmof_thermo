import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import gzip
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
df["pore_diameters"] = df["pore_diameters"].str[1]

counts = df["synthesizable"].value_counts()
print(counts)
print(df)
counts = df["synthesizable"].value_counts()
print(counts)

# Separate data by synthesizable status
df_synth = df[df["synthesizable"] == True]
df_not_synth = df[df["synthesizable"] == False]

# 5) Create hexbin plot
fig, ax = plt.subplots(figsize=(9, 6))

# Create hexbin for not synthesized (plot first so synthesized appears on top)
if len(df_not_synth) > 0:
    hb_red = ax.hexbin(df_not_synth["pore_diameters"], df_not_synth["ehull"], 
                    gridsize=50, cmap='Reds', alpha=0.6, mincnt=1,
                    extent=[0, 46, 0, 0.82], norm=LogNorm())

# Create hexbin for synthesized
if len(df_synth) > 0:
    hb_blue = ax.hexbin(df_synth["pore_diameters"], df_synth["ehull"], 
                    gridsize=50, cmap='Blues', alpha=0.6, mincnt=1,
                    extent=[0, 46, 0, 0.82], norm=LogNorm())

vmin_val = min(hb_red.get_array().min(), hb_blue.get_array().min())
vmax_val = max(hb_red.get_array().max(), hb_blue.get_array().max())

if len(df_not_synth) > 0:
    cbar1 = plt.colorbar(hb_red, ax=ax, pad=-0.02)
#    cbar1.set_label('Not Synthesized', fontsize=12, color='darkred')
    cbar1.ax.tick_params(labelsize=20)
    cbar1.mappable.set_clim(vmin=vmin_val, vmax=vmax_val)

if len(df_synth) > 0:
    # Position the second colorbar to the right of the first
    cbar2 = plt.colorbar(hb_blue, ax=ax, pad=0.01)  # Increased pad to avoid overlap
#    cbar2.set_label('Synthesized', fontsize=12, color='darkblue')
    cbar2.ax.tick_params(labelsize=20)
    cbar2.mappable.set_clim(vmin=vmin_val, vmax=vmax_val)


plt.xlabel("Pore Limiting Diameter (Å)", fontsize = 22)
ax.xaxis.set_label_coords(0.6, -0.1)
plt.ylabel("Δ$E_{\mathrm{hull}}$ (eV/atom)", fontsize = 22)



ax.tick_params(
    which='major', direction='in', length=26, width=1.8
)
ax.tick_params(
    which='minor', direction='in', length=14, width=1.6
)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(labelsize=20)
ax.minorticks_on()

# Create custom legend
plt.legend(handles=[
    plt.Line2D([0], [0], marker='h', color='w', markerfacecolor='#6495ED', markersize=12, label='Synthesized'),
    plt.Line2D([0], [0], marker='h', color='w', markerfacecolor='#FF6347', markersize=12, label='Hypothetical')
],
fontsize=20,
loc='upper right',
frameon=False
)

plt.xlim([0, 46])
plt.ylim([0, 0.82])

ax.text(
    -0.16,    # x position, just left of the left spine
    1.06,     # y position, just above the top spine
    "A",      # the label
    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
    fontsize=24,              # size of the letter
    fontweight="bold",        # make it bold
    va="top",                 # vertical alignment
    ha="left"                 # horizontal alignment
)

plt.tight_layout()
plt.savefig("FigureS16a")
plt.show()

