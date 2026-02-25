import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gzip

# Load data once
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

print(len(data))

# Build DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)

# Create figure with 2 subplots side by side
# Combined width = 2.8 + 2.8 = 5.6, height = max(2.5, 2.55) = 2.55
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.6, 2.55))

# ============== PLOT A (Formation Energy) ==============
df_plot1 = df[["formation_energy", "synthesizable", "chemsys", "frac_composition"]].copy()
df_plot1["group"] = "All MOFs"

sns.violinplot(
    data=df_plot1,
    x="group",
    y="formation_energy",
    hue="synthesizable",
    hue_order=[True, False],
    split=True,
    inner="quart",
    palette={True: "#6495ED", False: "#FF6347"},
    ax=ax1
)

ax1.set_xticklabels([])
ax1.tick_params(which='major', direction='in', length=10, width=1.25)
ax1.tick_params(which='minor', direction='in', length=5, width=1.25)
ax1.tick_params(axis='y', labelsize=9)
ax1.tick_params(axis='x', which='both', length=0)
ax1.minorticks_on()
ax1.set_ylim([-1.5, 0.65])
ax1.set_xlabel("", fontsize=1)
ax1.set_ylabel("Δ$E_{\mathrm{form}}$ (eV/atom)", fontsize=10)

for spine in ax1.spines.values():
    spine.set_linewidth(1.25)

ax1.text(-0.265, 1.00, "A", transform=ax1.transAxes, fontsize=11, 
         fontweight="bold", va="top", ha="left")
ax1.get_legend().remove()

# ============== PLOT B (E_hull) ==============
df_plot2 = df[["ehull", "synthesizable", "chemsys", "frac_composition"]].copy()
df_plot2["group"] = "All MOFs"

print(df_plot2["synthesizable"].value_counts())

sns.violinplot(
    data=df_plot2,
    x="group",
    y="ehull",
    hue="synthesizable",
    hue_order=[True, False],
    split=True,
    inner="quart",
    palette={False: "#FF6347", True: "#6495ED"},
    ax=ax2
)

ax2.set_xticklabels([])
ax2.tick_params(which='major', direction='in', length=10, width=1.25)
ax2.tick_params(which='minor', direction='in', length=5, width=1.25)
ax2.tick_params(axis='y', labelsize=9)
ax2.tick_params(axis='x', which='both', length=0)
ax2.minorticks_on()
ax2.set_ylim([0, 0.8])
ax2.set_xlabel("", fontsize=22)
ax2.set_ylabel("$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=10)

for spine in ax2.spines.values():
    spine.set_linewidth(1.25)

ax2.legend(title="Synthesized", fontsize=9, title_fontsize=9, 
           frameon=False, loc="upper right")

ax2.text(-0.22, 1.00, "B", transform=ax2.transAxes, fontsize=11,
         fontweight="bold", va="top", ha="left")

plt.tight_layout()
plt.savefig("Figure1.png", dpi=1500, bbox_inches='tight')
plt.show()
