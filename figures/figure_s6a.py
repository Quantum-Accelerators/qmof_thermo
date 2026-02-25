import json
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import gzip

temp = 298

# --- 1. Load the JSON data ---
# Replace 'data.json' with the path to your file
with gzip.open('UMA_MOF_phonon_results_with_gas_'+str(temp)+'.json.gz', 'rt') as f:
    data = json.load(f)

# --- 2. Extract all ehull_difference values ---
Ghull = [
    entry.get('Ghull')
    for entry in data.values()
]

Ehull = [
    entry.get('Ehull')
    for entry in data.values()
]

print(len(Ghull))

fig, ax = plt.subplots(figsize=(8, 7))
plt.scatter(Ehull, Ghull, s=50, c='k', alpha=0.7)


plt.plot([0.1, 0.6], [0.1, 0.6], linestyle='--', color='black', linewidth=2)

# Labels and styling
plt.xlabel(r'Δ$E_{\mathrm{hull}}$ (eV/atom)', fontsize=22)
plt.ylabel(r'Δ$G_{\mathrm{hull}}$ (eV/atom)', fontsize=22)

ax.minorticks_on()

ax.tick_params(
    which='major', direction='in', length=26, width=1.8
)
ax.tick_params(
    which='minor', direction='in', length=14, width=1.6
)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(labelsize=20)

#ax.minorticks_on()

ax.tick_params(
    axis='x',        # apply to the x axis
   # which='both',    # both major and minor ticks
   # length=0,
    pad = 10   # tick mark length = 0
)

plt.grid(False)

plt.xlim([0.1, 0.525])
plt.ylim([0.1, 0.525])

ticks = np.arange(0.1, 0.6, 0.1)  # or 0.1 step: np.arange(0.1, 0.6, 0.1)
ax.set_xticks(ticks)
ax.set_yticks(ticks)

plt.legend(title = str(temp)+' K', title_fontsize = 30, frameon=False, loc = "center right")

ax.text(
    -0.175,    # x position, just left of the left spine
    1.05,     # y position, just above the top spine
    "A",      # the label
    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
    fontsize=24,              # size of the letter
    fontweight="bold",        # make it bold
    va="top",                 # vertical alignment
    ha="left"                 # horizontal alignment
)


plt.tight_layout()
plt.savefig("FigureS6a")
plt.show()
# plt.savefig('ehull_difference_histogram.png', dpi=300)

