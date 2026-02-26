import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import gzip
# 1) Load your results JSON
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

with gzip.open("All_AI_MOF_results.json.gz", "rt") as f:
    data2 = json.load(f)

for qmof_id, properties in data2.items():
    if "model" in properties and "source" not in properties:
        properties["source"] = properties["model"]


data.update(data2)
#GHP-MOFassemble    
# 2) Build a tidy DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)

rename_map = {
    "Haranczyk_MOF5": "Mail-Order\nMOF-5",
    "Haranczyk_MOF74":  "hMOF-74",
    "ghp": "GHP-\nMOF-\nassemble",
    "Anderson": "Anderson-\nZr/Hf",
    "ToBaCCo": "ToBaCCo-\nCu",
    "MOFFUSION": "MOF-\nFUSION",
    "BoydWoo": "Boyd-\nWoo"
    # add as many as you like…
}

# 3) apply it
df["source"] = df["source"].replace(rename_map)

df = df[["ehull", "source", "synthesizable", "chemsys"]]
#df = df[df["chemsys"] == "C-H-N-O-Zn"]
# 3) Choose which sources you want to plot (e.g. CSD, CoRE, BoydWoo)
filter_list = ["Pyrene", "CSD", "CoRE", "hMOF-74", "Mail-Order\nMOF-5", "GMOF", "Boyd-\nWoo", "Anderson-\nZr/Hf", "ToBaCCo-\nCu", "GHP-\nMOF-\nassemble", "MOF-\nFUSION", "MOFDiff"]

# 4) Filter the DataFrame to only those sources

counts = df["synthesizable"].value_counts()
print(counts)
counts = df["source"].value_counts()
print(counts)


df = df[df["source"].isin(filter_list)].copy()
df["source"] = pd.Categorical(
    df["source"],
    categories=filter_list,  # exact same strings & order
    ordered=True
)

blue_shades = sns.color_palette("Blues", 3)
red_shades  = sns.color_palette("Reds", 6)
green_shades = sns.color_palette("Greens", 3)
palette = dict(zip(filter_list, blue_shades + red_shades + green_shades))


fig, ax = plt.subplots(figsize=(6.75, 3))

sns.violinplot(
    data=df,
    x="source",
    hue="source",
    y="ehull",
#    hue="synthesizable",
#    split=True,
#    hue_order = [True, False],
    inner="quart",
    order = filter_list,
    ax =ax,
    width = 0.9,
    palette = palette,
    legend=False,
 #    palette={ False: "#FF6347", True: "#6495ED"},
    cut=0,
    linewidth = 0.9
)

ax.set_xticks(np.arange(len(filter_list)))

labels = ax.get_xticklabels()
for lbl in labels:
    lbl.set_fontsize(9)
ax.set_xticklabels(labels)


xticks = ax.get_xticks()                # array([0., 1., 2., …])
xlabels = [lbl.get_text() for lbl in ax.get_xticklabels()]



i_pyrene  = xlabels.index("CoRE")
i_boydwoo = xlabels.index("hMOF-74")    # match exact casing!

# 5) Compute the halfway‐point and draw the dashed line
x_line = (xticks[i_pyrene] + xticks[i_boydwoo]) / 2
ax.axvline(x=x_line, linestyle="--", color="black", linewidth=1.25)


i_pyrene  = xlabels.index("ToBaCCo-\nCu")
i_boydwoo = xlabels.index("GHP-\nMOF-\nassemble")    # match exact casing!

# 5) Compute the halfway‐point and draw the dashed line
x_line = (xticks[i_pyrene] + xticks[i_boydwoo]) / 2
ax.axvline(x=x_line, linestyle="--", color="black", linewidth=1.25)



new_labels = [
    f"{src} ({counts[src]})" 
    if src in ["ToBaCCo-\nCu"]#, "Anderson-\nZr/Hf"]#, "GHP-\nMOFassemble"]
    else f"{src}\n({counts[src]})"
    for src in filter_list
]

ax.set_xticklabels(new_labels, rotation=22.5, fontsize=8)

labels = ax.get_xticklabels()
ghp_index = filter_list.index("GHP-\nMOF-\nassemble")
labels[ghp_index].set_verticalalignment('top')
labels[ghp_index].set_y(0.03)

labels = ax.get_xticklabels()
ghp_index = filter_list.index("Mail-Order\nMOF-5")
labels[ghp_index].set_verticalalignment('top')
labels[ghp_index].set_y(0.015)

labels = ax.get_xticklabels()
ghp_index = filter_list.index("Boyd-\nWoo")
labels[ghp_index].set_verticalalignment('top')
labels[ghp_index].set_y(-0.015)

labels = ax.get_xticklabels()
ghp_index = filter_list.index("Anderson-\nZr/Hf")
labels[ghp_index].set_verticalalignment('top')
labels[ghp_index].set_y(0.02)

labels = ax.get_xticklabels()
ghp_index = filter_list.index("ToBaCCo-\nCu")
labels[ghp_index].set_verticalalignment('top')
labels[ghp_index].set_y(0.0)






ax.tick_params(
    which='major', direction='in', length=15, width=1.5
)
ax.tick_params(
    which='minor', direction='in', length=14, width=1.6
)

for spine in ax.spines.values():
    spine.set_linewidth(1.25)

ax.tick_params(
    axis='y',        # apply to the x axis
    which='both',    # both major and minor ticks
    length=0         # tick mark length = 0
)

ax.tick_params(axis='y', labelsize=10)

plt.ylabel("$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize = 10)
plt.xlabel("")

plt.ylim([-.10, 0.8])

for y in ax.get_yticks():
    ax.axhline(y=y, linestyle="-", linewidth=1.05, color="grey", alpha=0.8, zorder=0)



plt.tight_layout()
plt.savefig("Figure4", dpi=1500, bbox_inches='tight')
plt.show()
