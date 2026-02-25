import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import norm
from scipy.stats import shapiro
from scipy.stats import normaltest
from scipy.stats import anderson
from scipy.stats import kstest
import gzip
#ehull = 0.2

with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# 2) Build a tidy DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)

df = df[[ "ehull", "synthesizable"]]
df = df[df["synthesizable"] == True]

counts = df["synthesizable"].value_counts()
print(counts)

ehulls = df["ehull"].values

fig, ax = plt.subplots(figsize=(8, 6))
stats.probplot(ehulls, dist="norm", plot=plt)


dots_line, fit_line = ax.get_lines()

# change dot size & color
dots_line.set_marker('o')             # ensure marker is a circle
dots_line.set_markersize(5)           # e.g. size 8
dots_line.set_markeredgewidth(1.5)    # edge width
dots_line.set_markeredgecolor('k')    # black edge
dots_line.set_markerfacecolor('k')   # e.g. matplotlib “C1” (orange) fill
dots_line.set_label('MOF Quantiles')
# change fit‐line width & color (and style, if you like)
fit_line.set_linewidth(2.5)
fit_line.set_color('r')              # e.g. “C3” (green)
fit_line.set_linestyle('-')          # dashed
fit_line.set_label('Normal Quantiles')

plt.title("")

plt.legend(fontsize = 22,
loc = 'best',
frameon = False)

ax.tick_params(
    which='major', direction='in', length=26, width=1.8
)
ax.tick_params(
    which='minor', direction='in', length=14, width=1.6
)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(labelsize=20)

ax.minorticks_on()

#plt.ylim([0, 0.55])
#plt.title("Q–Q Plot of Synthesized Ehull vs. Normal Distribution")
plt.xlabel("Normal Theoretical Quantiles", fontsize = 22)
plt.ylabel("Δ$E_{\mathrm{hull}}$ (eV/atom)", fontsize = 22)

#ax.text(
#    -0.14,    # x position, just left of the left spine
#    1.06,     # y position, just above the top spine
#    "B",      # the label
#    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
#    fontsize=24,              # size of the letter
#    fontweight="bold",        # make it bold
#    va="top",                 # vertical alignment
#    ha="left"                 # horizontal alignment
#)

plt.tight_layout()
plt.savefig("FigureS10")
plt.show()


