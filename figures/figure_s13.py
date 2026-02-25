import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import gzip

# 1) Load your results JSON (no AI MOFs)
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# 2) Build a tidy DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)

rename_map = {
    "Haranczyk_MOF5": "Mail-Order\nMOF-5",
    "Haranczyk_MOF74":  "hMOF-74",
    "ghp": "GHP-\nMOFassemble",
}

# 3) Apply renaming
df["source"] = df["source"].replace(rename_map)

# Extract first pore diameter
def get_pore_diameter(pore_list):
    if isinstance(pore_list, list) and len(pore_list) > 0:
        return pore_list[0]
    return None

df["pore_diameter"] = df["pore_diameters"].apply(get_pore_diameter)

# Keep only entries with pore diameter data
df = df[df["pore_diameter"].notna()].copy()

print(f"Entries with pore diameter data: {len(df)}")

df = df[["pore_diameter", "source", "synthesizable", "chemsys"]]

# 3) Choose which sources you want to plot
filter_list = ["Pyrene", "CSD", "CoRE", "hMOF-74", "Mail-Order\nMOF-5", "GMOF", 
               "BoydWoo", "Anderson", "ToBaCCo"]  # Removed AI MOFs

# 4) Filter the DataFrame to only those sources
counts = df["synthesizable"].value_counts()
print(counts)
counts = df["source"].value_counts()
print(counts)

df = df[df["source"].isin(filter_list)].copy()
df["source"] = pd.Categorical(
    df["source"],
    categories=filter_list,
    ordered=True
)

blue_shades = sns.color_palette("Blues", 3)
red_shades = sns.color_palette("Reds", 6)
palette = dict(zip(filter_list, blue_shades + red_shades))

fig, ax = plt.subplots(figsize=(6.75, 3))

sns.violinplot(
    data=df,
    x="source",
    hue="source",
    y="pore_diameter",  # Changed
    inner="quart",
    order=filter_list,
    ax=ax,
    width=0.9,
    palette=palette,
    legend=False,
    cut=0,
    linewidth=0.9
)

ax.set_xticks(np.arange(len(filter_list)))

labels = ax.get_xticklabels()
for lbl in labels:
    lbl.set_fontsize(9)
ax.set_xticklabels(labels)

xticks = ax.get_xticks()
xlabels = [lbl.get_text() for lbl in ax.get_xticklabels()]

i_pyrene = xlabels.index("CoRE")
i_boydwoo = xlabels.index("hMOF-74")

x_line = (xticks[i_pyrene] + xticks[i_boydwoo]) / 2
ax.axvline(x=x_line, linestyle="--", color="black", linewidth=1.25)

new_labels = [
    f"{src} ({counts[src]})"
    for src in filter_list
]
ax.set_xticklabels(new_labels, rotation=22.5, fontsize=8)

labels = ax.get_xticklabels()
mof5_index = filter_list.index("Mail-Order\nMOF-5")
labels[mof5_index].set_verticalalignment('top')
labels[mof5_index].set_y(0.015)

ax.tick_params(which='major', direction='in', length=15, width=1.5)
ax.tick_params(which='minor', direction='in', length=14, width=1.6)

for spine in ax.spines.values():
    spine.set_linewidth(1.25)

ax.tick_params(axis='y', which='both', length=0)
ax.tick_params(axis='y', labelsize=10)

plt.ylabel("Largest Cavity Diameter (Å)", fontsize=10)  # Changed label
plt.xlabel("")

plt.ylim([-5, 50])  # Typical pore diameter range

for y in ax.get_yticks():
    ax.axhline(y=y, linestyle="-", linewidth=1.05, color="grey", alpha=0.8, zorder=0)

plt.tight_layout()
plt.savefig("FigureS13", dpi=1500, bbox_inches='tight')
plt.show()
