import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import gzip

with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# Create 3x2 subplot grid
fig, axes = plt.subplots(3, 2, figsize=(14, 12))

elements = ["H", "O", "N"]
labels = ["A", "B", "C"]
element_names = {
    "H": "Hydrogen",
    "O": "Oxygen", 
    "N": "Nitrogen"
}

# Build base dataframe once
df_all = pd.DataFrame.from_dict(data, orient="index").rename_axis("qmof_id")

for row, (mol, label) in enumerate(zip(elements, labels)):
    # Filter for this element
    df = df_all[df_all["frac_composition"].apply(lambda d: mol in d)]
    df = df.assign(molfrac=df["frac_composition"].apply(lambda d: d[mol]))
    
    # Left column: formation_energy
    ax_left = axes[row, 0]
    hb_left = ax_left.hexbin(df["molfrac"], df["formation_energy"],
                             gridsize=50, cmap='viridis', mincnt=1,
                             norm=LogNorm())
    
    # Colorbar for left
    cb_left = plt.colorbar(hb_left, ax=ax_left, pad=0.02)
    cb_left.ax.tick_params(labelsize=16)
    
    # Styling for left
    ax_left.tick_params(which='major', direction='in', length=20, width=1.5)
    ax_left.tick_params(which='minor', direction='in', length=10, width=1.2)
    ax_left.tick_params(labelsize=16)
    ax_left.minorticks_on()
    for spine in ax_left.spines.values():
        spine.set_linewidth(2)
    
    # Y-label only on left column
    ax_left.set_ylabel("Δ$E_{\mathrm{form}}$ (eV/atom)", fontsize=18)
    
    # X-label only on bottom row
    if row == 2:
        ax_left.set_xlabel(f"{element_names[mol]} Fraction", fontsize=18)
    
    # Add A, B, C labels to left column only
    ax_left.text(-0.165, 1.06, label, transform=ax_left.transAxes,
                fontsize=20, fontweight="bold", va="top", ha="left")
    
    # Right column: ehull
    ax_right = axes[row, 1]
    hb_right = ax_right.hexbin(df["molfrac"], df["ehull"],
                               gridsize=50, cmap='viridis', mincnt=1,
                               norm=LogNorm())
    
    # Colorbar for right
    cb_right = plt.colorbar(hb_right, ax=ax_right, pad=0.02)
    cb_right.ax.tick_params(labelsize=16)
    
    # Styling for right
    ax_right.tick_params(which='major', direction='in', length=20, width=1.5)
    ax_right.tick_params(which='minor', direction='in', length=10, width=1.2)
    ax_right.tick_params(labelsize=16)
    ax_right.minorticks_on()
    for spine in ax_right.spines.values():
        spine.set_linewidth(2)
    
    # Y-label for right column
    ax_right.set_ylabel("Δ$E_{\mathrm{hull}}$ (eV/atom)", fontsize=18)
    
    # X-label only on bottom row
    if row == 2:
        ax_right.set_xlabel(f"{element_names[mol]} Fraction", fontsize=18)

plt.tight_layout()
plt.savefig("FigureS24.png", bbox_inches='tight', dpi=500)
plt.show()
