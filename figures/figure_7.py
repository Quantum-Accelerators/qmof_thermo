from __future__ import annotations

import gzip
import json
import re
from collections import Counter

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from ase import Atoms
from matplotlib.colors import LogNorm
from scipy.stats import spearmanr

# Load data
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# ==================== DATA PREPARATION FOR PLOT A ====================
counter = Counter()
for entry in data.values():
    decomp = entry.get("decomposition_products", {})
    counter.update(decomp.keys())

top10 = counter.most_common(10)
products, counts = zip(*top10, strict=False)

df_a = pd.DataFrame({"product": products, "count": counts})


def to_hill_subscript(prod_str):
    custom_formulas = {
        "H3N": "NH3",
        "ZnO": "ZnO",
        "H4N2": "N2H4",
        "H6C(NO)2": "[NH4+][NH2COO–]",
        "H5NO": "NH3·H2O",
        "H5CNO3": "[NH4+][HCO3–]",
    }

    if prod_str in custom_formulas:
        formula = custom_formulas[prod_str]
    else:
        formula = Atoms(prod_str).get_chemical_formula()

    formula = re.sub(r"(\d+)", r"$_{\1}$", formula)
    return re.sub(r"([+\-–])", r"$^{\1}$", formula)


df_a["pretty"] = df_a["product"].map(to_hill_subscript)

# ==================== DATA PREPARATION FOR PLOTS B & C ====================
mol = "C"
df = pd.DataFrame.from_dict(data, orient="index").rename_axis("qmof_id")
df = df[df["frac_composition"].apply(lambda d: mol in d)]
df = df.assign(molfrac=df["frac_composition"].apply(lambda d: d[mol]))

corr_eform, pval_eform = spearmanr(df["molfrac"], df["formation_energy"])
corr_ehull, pval_ehull = spearmanr(df["molfrac"], df["ehull"])

print("\n=== Spearman Correlation Coefficients ===")
print(
    f"Carbon Fraction vs Formation Energy: ρ = {corr_eform:.4f}, p-value = {pval_eform}"
)
print(
    f"Carbon Fraction vs Energy Above Hull: ρ = {corr_ehull:.4f}, p-value = {pval_ehull}"
)
print(f"Sample size: {len(df)}")

# Get shared colorbar range
fig_temp, (ax1_temp, ax2_temp) = plt.subplots(1, 2)
hb1_temp = ax1_temp.hexbin(
    df["molfrac"],
    df["formation_energy"],
    gridsize=80,
    mincnt=1,
    extent=[0, 0.7, -2.3, 0.7],
)
hb2_temp = ax2_temp.hexbin(
    df["molfrac"], df["ehull"], gridsize=80, mincnt=1, extent=[0, 0.7, -0.1, 0.85]
)
vmin = min(hb1_temp.get_array().min(), hb2_temp.get_array().min())
vmax = max(hb1_temp.get_array().max(), hb2_temp.get_array().max())
plt.close(fig_temp)

# ==================== CREATE COMBINED FIGURE ====================
# Use 4 columns so top plot can be centered and narrower
fig = plt.figure(figsize=(7.45, 5.5))
gs = gridspec.GridSpec(
    2, 14, figure=fig, height_ratios=[1, 1], wspace=0.45, hspace=0.25
)

# ============== PLOT A: Top Decomposition Products (centered, spans 2 middle columns) ==============
ax1 = fig.add_subplot(gs[0, 5:11])  # Columns 1 and 2 (middle)

sns.barplot(data=df_a, x="count", y="product", palette="viridis", ax=ax1)

ax1.tick_params(labelsize=8)
ax1.minorticks_on()
ax1.tick_params(which="major", direction="in", length=10, width=1.25)
ax1.tick_params(which="minor", direction="in", length=5, width=1.25)
ax1.tick_params(axis="y", which="both", length=0)

for spine in ax1.spines.values():
    spine.set_linewidth(1.25)

ax1.set_xlim([0, 21000])
ax1.set_xlabel("Number of Occurrences", fontsize=8)
ax1.set_ylabel("Decomposition Product", fontsize=8)
ax1.set_yticklabels(df_a["pretty"], fontsize=8)

ax1.text(
    -0.5,
    1,
    "A",
    transform=ax1.transAxes,
    fontsize=11,
    fontweight="bold",
    va="top",
    ha="left",
)

# ============== PLOT B: Formation Energy vs Carbon Fraction ==============
ax2 = fig.add_subplot(gs[1, 0:6])  # Columns 0 and 1

hb2 = ax2.hexbin(
    df["molfrac"],
    df["formation_energy"],
    gridsize=80,
    cmap="viridis",
    mincnt=1,
    extent=[0, 0.7, -2.3, 0.7],
    norm=LogNorm(vmin=vmin, vmax=vmax),
)

ax2.text(
    0.95,
    0.2,
    f"$ρ$ = {corr_eform:.3f}",
    transform=ax2.transAxes,
    fontsize=8,
    verticalalignment="bottom",
    horizontalalignment="right",
    bbox={"boxstyle": "round", "facecolor": "white", "alpha": 1},
)

ax2.tick_params(which="major", direction="in", length=10, width=1.25)
ax2.tick_params(which="minor", direction="in", length=5, width=1.25)
ax2.tick_params(labelsize=8)
ax2.minorticks_on()

for spine in ax2.spines.values():
    spine.set_linewidth(1.25)

ax2.set_xlim([0, 0.675])
ax2.set_ylim([-2.3, 0.7])
ax2.set_xlabel("Carbon Fraction in MOF", fontsize=8)
ax2.set_ylabel(r"$ΔE_{\mathrm{form}}$ (eV/atom)", fontsize=8)

ax2.text(
    -0.2,
    1.00,
    "B",
    transform=ax2.transAxes,
    fontsize=11,
    fontweight="bold",
    va="top",
    ha="left",
)

# ============== PLOT C: Ehull vs Carbon Fraction ==============
ax3 = fig.add_subplot(gs[1, 7:14])  # Columns 2 and 3

hb3 = ax3.hexbin(
    df["molfrac"],
    df["ehull"],
    gridsize=80,
    cmap="viridis",
    mincnt=1,
    extent=[0, 0.7, -0.1, 0.85],
    norm=LogNorm(vmin=vmin, vmax=vmax),
)

ax3.text(
    0.95,
    0.95,
    f"$ρ$ = {corr_ehull:.3f}",
    transform=ax3.transAxes,
    fontsize=8,
    verticalalignment="top",
    horizontalalignment="right",
    bbox={"boxstyle": "round", "facecolor": "white", "alpha": 1},
)

cb = plt.colorbar(hb3, ax=ax3, pad=0.02)
cb.ax.tick_params(labelsize=8)

ax3.tick_params(which="major", direction="in", length=10, width=1.25)
ax3.tick_params(which="minor", direction="in", length=5, width=1.25)
ax3.tick_params(labelsize=8)
ax3.minorticks_on()

for spine in ax3.spines.values():
    spine.set_linewidth(1.25)

ax3.set_xlim([0, 0.675])
ax3.set_ylim([-0.05, 0.85])
ax3.set_xlabel("Carbon Fraction in MOF", fontsize=8)
ax3.set_ylabel(r"$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=8)

ax3.text(
    -0.17,
    1.0,
    "C",
    transform=ax3.transAxes,
    fontsize=11,
    fontweight="bold",
    va="top",
    ha="left",
)

plt.tight_layout()
plt.savefig("Figure7.png", dpi=1500, bbox_inches="tight")
plt.show()
