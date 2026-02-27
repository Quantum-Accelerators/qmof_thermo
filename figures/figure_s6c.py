from __future__ import annotations

import gzip
import json

import matplotlib.pyplot as plt
import numpy as np

temp = 298

# --- 1. Load the JSON data ---
# Replace 'data.json' with the path to your file
with gzip.open("UMA_MOF_phonon_results_with_gas_" + str(temp) + ".json.gz", "rt") as f:
    data = json.load(f)

# --- 2. Extract all ehull_difference values ---
Gform = [entry.get("Gform") for entry in data.values()]

Eform = [entry.get("Eform") for entry in data.values()]

with gzip.open("UMA_MP_phonon_results_with_gas_" + str(temp) + ".json.gz", "rt") as f:
    data = json.load(f)

# --- 2. Extract all ehull_difference values ---
mp_Gform = [
    entry.get("Gform")
    for entry in data.values()
    if entry.get("Ghull") == 0 or entry.get("Ehull") == 0
]

mp_Eform = [
    entry.get("Eform")
    for entry in data.values()
    if entry.get("Ehull") == 0 or entry.get("Ghull") == 0
]


fig, ax = plt.subplots(figsize=(8, 7))
plt.scatter(Eform, Gform, s=50, c="k", alpha=0.7, label="MOFs")

plt.scatter(mp_Eform, mp_Gform, s=50, c="r", alpha=0.7, label="Hull Structures")

plt.plot([-2, 0.75], [-2, 0.75], linestyle="--", color="black", linewidth=2)

# Labels and styling
plt.xlabel(r"Δ$E_{\mathrm{form}}$ (eV/atom)", fontsize=22)
plt.ylabel(r"Δ$G_{\mathrm{form}}$ (eV/atom)", fontsize=22)

ax.minorticks_on()

ax.tick_params(which="major", direction="in", length=26, width=1.8)
ax.tick_params(which="minor", direction="in", length=14, width=1.6)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(labelsize=20)

# ax.minorticks_on()

ax.tick_params(
    axis="x",  # apply to the x axis
    # which='both',    # both major and minor ticks
    # length=0,
    pad=10,  # tick mark length = 0
)

plt.grid(False)

plt.xlim([-1.5, 0.25])
plt.ylim([-1.5, 0.25])

ticks = np.arange(-1.5, 0.25, 0.2)  # or 0.1 step: np.arange(0.1, 0.6, 0.1)
ax.set_xticks(ticks)
ax.set_yticks(ticks)


legend1 = plt.legend(
    fontsize=18,
    frameon=True,
    loc="lower right",
    bbox_to_anchor=(1, 0.1),
    edgecolor="black",
)

# Add second legend for temperature (as text box)
ax.add_artist(legend1)  # Keep first legend
ax.text(
    0.88,
    0.5,
    str(temp) + " K",
    transform=ax.transAxes,
    fontsize=30,
    verticalalignment="center",
    horizontalalignment="center",
)  # ,
#  bbox=dict(boxstyle='round', facecolor='white', alpha=1, edgecolor='black'))


ax.text(
    -0.175,  # x position, just left of the left spine
    1.05,  # y position, just above the top spine
    "C",  # the label
    transform=ax.transAxes,  # interpret x,y in axis fraction (0 to 1)
    fontsize=24,  # size of the letter
    fontweight="bold",  # make it bold
    va="top",  # vertical alignment
    ha="left",  # horizontal alignment
)

plt.tight_layout()
plt.savefig("FigureS6c")
plt.show()
# plt.savefig('ehull_difference_histogram.png', dpi=300)
