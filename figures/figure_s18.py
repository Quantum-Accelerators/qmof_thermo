from __future__ import annotations

import gzip
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)

# 1) Load your results JSON
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# 2) Build a tidy DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
    .reset_index()
    .rename(columns={"index": "qmof_id"})
)


# Extract first pore diameter
def get_first_pore_diameter(pore_diameters):
    """Extract the first pore diameter from the list"""
    if isinstance(pore_diameters, list) and len(pore_diameters) > 0:
        return pore_diameters[0]
    return np.nan


# Add pore diameter column
df["pore_diameter"] = df["pore_diameters"].apply(get_first_pore_diameter)

# Keep relevant columns
df = df[["ehull", "topology", "synthesizable", "chemsys", "pore_diameter"]]

# 3) Choose which topologies you want to plot
filter_list = [
    "sql",
    "pcu",
    "hcb",
    "dia",
    "kgd",
    "fes",
    "rtl",
    "bey",
    "pts",
    "bex",
    "fsc",
    "rna",
    "cds",
    "dmd",
    "ant",
]

# 4) Filter the DataFrame
df = df[df["synthesizable"]]
counts = df["synthesizable"].value_counts()
print(counts)
counts = df["topology"].value_counts()
print(counts)

median_ehull = df.groupby("topology")["ehull"].median().sort_values()
counts = df.groupby("topology").size()
median_pore = df.groupby("topology")["pore_diameter"].median().sort_values()
sorted_filter_list_ehull = median_ehull.index.tolist()

print("All topology median ehull:")
for topology in sorted_filter_list_ehull:
    median_val = median_ehull[topology]
    count = counts[topology]
    print(f"  {topology}: {median_val:.4f} {count}")


# Create a combined statistics DataFrame
stats_df = pd.DataFrame(
    {
        "topology": median_ehull.index,
        "median_ehull": median_ehull.to_numpy(),
        "count": [counts[top] for top in median_ehull.index],
        "median_pore_diameter": [
            median_pore[top] if top in median_pore.index else np.nan
            for top in median_ehull.index
        ],
    }
)

# Sort by median ehull (already sorted from median_ehull)
stats_df = stats_df.reset_index(drop=True)

# Export to Excel
stats_df.to_excel("topology_statistics.xlsx", index=False)
print("\nStatistics exported to 'topology_statistics.xlsx'")


df = df[df["topology"].isin(filter_list)].copy()

# 5) AUTO-SORT by median ehull
median_ehull = df.groupby("topology")["ehull"].median().sort_values()
sorted_filter_list_ehull = median_ehull.index.tolist()

print("Topologies sorted by median ehull:")
for topology in sorted_filter_list_ehull:
    median_val = median_ehull[topology]
    print(f"  {topology}: {median_val:.4f}")

# 6) AUTO-SORT by median pore diameter
median_pore_diameter = df.groupby("topology")["pore_diameter"].median().sort_values()
sorted_filter_list_pore = median_pore_diameter.index.tolist()

print("\nTopologies sorted by median pore diameter:")
for topology in sorted_filter_list_pore:
    median_val = median_pore_diameter[topology]
    print(f"  {topology}: {median_val:.4f}")

# 7) Use ehull sorting for the plot (you can change this to sorted_filter_list_pore if preferred)
sorted_filter_list = sorted_filter_list_ehull

# Set categorical order using sorted list
df["topology"] = pd.Categorical(
    df["topology"],
    categories=sorted_filter_list,  # Use sorted order
    ordered=True,
)

fig, ax = plt.subplots(figsize=(14, 6))
sns.violinplot(
    data=df,
    x="topology",
    hue="topology",
    y="ehull",
    inner="quart",
    order=sorted_filter_list,  # Use sorted order
    ax=ax,
    width=0.9,
    legend=False,
    cut=0,
)

ax.set_xticks(np.arange(len(sorted_filter_list)))
labels = ax.get_xticklabels()
for lbl in labels:
    lbl.set_fontsize(22)
#  lbl.set_fontweight('bold')
ax.set_xticklabels(labels)

# Update labels with counts using sorted order
new_labels = [f"{src} ({counts[src]})" for src in sorted_filter_list]
ax.set_xticklabels(new_labels, rotation=20, fontsize=18)

ax.tick_params(which="major", direction="in", length=26, width=1.8)
ax.tick_params(which="minor", direction="in", length=14, width=1.6)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(axis="y", which="both", length=0)

ax.tick_params(axis="y", labelsize=20)
plt.ylabel(r"Δ$E_{\mathrm{hull}}$ (eV/atom)", fontsize=22)
plt.xlabel("")
plt.ylim([-0.10, 0.8])

for y in ax.get_yticks():
    ax.axhline(y=y, linestyle="-", linewidth=1.3, color="grey", zorder=0)

plt.tight_layout()
plt.savefig("FigureS18", dpi=1000)

plt.show()
