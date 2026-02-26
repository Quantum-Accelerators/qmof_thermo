import json
import pandas as pd
import matplotlib.pyplot as plt
from ase import Atoms
import re
import gzip

with gzip.open("MP_All_results.json.gz", "rt") as f:
    data = json.load(f)

mol = "C"

# build base df and extract co2_frac
df = (
    pd.DataFrame.from_dict(data, orient="index")
      .rename_axis("qmof_id")
      .reset_index()
)
df = df[df["frac_composition"].apply(lambda d: mol in d)]
df = df.assign(molfrac=df["frac_composition"].apply(lambda d: d[mol]))
df = df[["formula", "ehull", "formation_energy_per_atom", "molfrac"]]
df = df[df["ehull"] == 0.0000]

df['formula_hill'] = df['formula'].map(lambda f: Atoms(f).get_chemical_formula())

def latex_subscript(formula: str) -> str:
    # wrap each run of digits in $_{…}$
    return re.sub(r'(\d+)', r'$_{\1}$', formula)

df['formula_latex'] = df['formula_hill'].map(latex_subscript)

# specify which formulas to highlight and their colors
to_highlight = {
    "C":   "#FF4136",   # red
    "CO2": "#0074D9",   # blue
    "H4C": "#2ECC40",   # green  (or use "CH4" if your dataframe uses that spelling)
}

print(len(df))

fig, ax = plt.subplots(figsize=(3.4, 2.5))

# overplot each special formula
for formula, color in to_highlight.items():
    sub = df[df["formula"] == formula]
    if not sub.empty:
        ax.scatter(sub["molfrac"],
                   sub["formation_energy_per_atom"],
                   c=color, s=30, linewidth=1.2,
                   label = sub["formula_latex"].iloc[0])

rest = df[~df["formula"].isin(to_highlight)]
ax.scatter(rest["molfrac"],
           rest["formation_energy_per_atom"],
           c="black", s=8, alpha=0.7, edgecolors="none",
           label="Other")


# styling
ax.tick_params(which="major", direction="in", length=10, width=1.25)
ax.tick_params(which="minor", direction="in", length=5, width=1.25)
ax.tick_params(labelsize=8)
ax.minorticks_on()
for spine in ax.spines.values():
    spine.set_linewidth(1.25)

ax.set_xlabel("Carbon Fraction", fontsize=8)
plt.title("Only Hull Materials", fontsize = 10)
ax.set_ylabel("$ΔE_{\mathrm{form}}$ (eV/atom)",       fontsize=8)
ax.set_xlim([0, 1.02])
ax.set_ylim([-3, 0.1])

# legend: one entry per scatter call, placed best

leg = ax.legend(frameon=False, fontsize=8, loc="center right")
# you can also tweak marker scale in the legend if you like:
for handle in leg.legend_handles:
    handle._sizes = [30]


ax.text(
    -0.195,    # x position, just left of the left spine
    0.99,     # y position, just above the top spine
    "B",      # the label
    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
    fontsize=11,              # size of the letter
    fontweight="bold",        # make it bold
    va="top",                 # vertical alignment
    ha="left"                 # horizontal alignment
)


plt.tight_layout()
plt.savefig("FigureS23B.png", bbox_inches="tight", dpi=1500)
plt.show()

