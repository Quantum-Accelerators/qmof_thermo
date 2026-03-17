from __future__ import annotations

import gzip
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MultipleLocator
from scipy.stats import norm

with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# 2) Build a tidy DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
    .reset_index()
    .rename(columns={"index": "qmof_id"})
)

df = df[["ehull", "synthesizable"]]
df = df[df["synthesizable"]]

ehulls = df["ehull"].to_numpy()


mu, sigma = ehulls.mean(), ehulls.std()

print(f"the mean is {mu} eV/atom")
bin_width = 0.02

# --- 3. Plot the histogram ---
fig, ax = plt.subplots(figsize=(3.15, 2.25))
sns.histplot(
    ehulls,
    stat="probability",
    binwidth=bin_width,
    edgecolor="black",
    label="MOF Distribution",
)
plt.xlabel(r"$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=8)
plt.ylabel("Relative Frequency", fontsize=8)
# plt.title('Distribution of Synthesized Ehull')
plt.grid(True, linestyle="--", alpha=0.5)


x = np.linspace(ehulls.min(), ehulls.max(), 200)
pdf_prob = norm.pdf(x, loc=mu, scale=sigma) * bin_width

ax.plot(x, pdf_prob, color="red", linewidth=1.5, label="Normal Distribution")


ax.tick_params(which="major", direction="in", length=10, width=1.25)
ax.tick_params(which="minor", direction="in", length=5, width=1.25)

for spine in ax.spines.values():
    spine.set_linewidth(1.25)

ax.tick_params(labelsize=8)

ax.minorticks_on()

ax.tick_params(
    axis="x",  # apply to the x axis
    which="both",  # both major and minor ticks
    length=0,  # tick mark length = 0
)

plt.grid(False)


# Add this after creating the plot and before plt.tight_layout()
ax.yaxis.set_major_locator(MultipleLocator(0.05))

ax.legend(fontsize=8, loc="best", frameon=False, bbox_to_anchor=(1.03, 1.02))

# ax.text(
#    -0.19,    # x position, just left of the left spine
#    1.0,     # y position, just above the top spine
#    "A",      # the label
#    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
#    fontsize=11,              # size of the letter
#    fontweight="bold",        # make it bold
#    va="top",                 # vertical alignment
#    ha="left"                 # horizontal alignment
# )
plt.ylim([0, 0.15])
plt.tight_layout()
plt.savefig("FigureS9", dpi=1500, bbox_inches="tight")
plt.show()
# plt.savefig('ehull_difference_histogram.png', dpi=300)
