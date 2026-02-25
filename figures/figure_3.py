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
import numpy as np
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

sorted_ehulls = np.sort(ehulls)

# Calculate percentiles (0 to 100)
percentiles = np.linspace(0, 1, len(sorted_ehulls))

# Create the plot
fig, ax = plt.subplots(figsize=(3.25, 2.5))
plt.plot(sorted_ehulls, percentiles, 'b-', linewidth=1.5)

#plt.grid(True, alpha=0.3)
percentile_markers = [
 #   (50, 'Median', 'green'),
    (0.80, '0.80', 'brown'),
    (0.90, '0.90', 'red'),
    (0.95, '0.95', 'purple'),
    (0.99, '0.99', 'green')
]

for percentile_val, label, color in percentile_markers:
    ehull_val = np.percentile(sorted_ehulls, 100*percentile_val)
    
    # Horizontal line from y-axis (x=0) to the curve (x=ehull_val)
    ax.plot([0, ehull_val], [percentile_val, percentile_val], 
            color=color, linestyle='--', linewidth=1, alpha=0.7)
    
    # Vertical line from the curve (y=percentile_val) down to x-axis (y=0)
    ax.plot([ehull_val, ehull_val], [0, percentile_val], 
            color=color, linestyle='--', linewidth=1, alpha=0.7,
            label=f'{label}')#: {ehull_val:.3f}')

# Add some useful reference lines

# Print some useful statistics
print(f"Number of MOFs: {len(sorted_ehulls)}")
print(f"Median ehull: {np.median(sorted_ehulls):.4f} eV/atom")
print(f"80th percentile: {np.percentile(sorted_ehulls, 80):.4f} eV/atom")
print(f"90th percentile: {np.percentile(sorted_ehulls, 90):.4f} eV/atom")
print(f"95th percentile: {np.percentile(sorted_ehulls, 95):.4f} eV/atom")
print(f"99th percentile: {np.percentile(sorted_ehulls, 99):.4f} eV/atom")
print(f"Max ehull: {np.max(sorted_ehulls):.4f} eV/atom")
plt.title("")

plt.legend(fontsize = 9,
#title = "Quantile",
loc = 'best',
frameon = False)

ax.tick_params(
    which='major', direction='in', length=15, width=1.25
)
ax.tick_params(
    which='minor', direction='in', length=8, width=1.25
)

for spine in ax.spines.values():
    spine.set_linewidth(1.25)

ax.tick_params(labelsize=9)

ax.minorticks_on()

ax.tick_params(axis="x", pad=4)

plt.ylim([-0.5/100, 100.5/100])
plt.xlim([0, 0.8])
#plt.title("Q–Q Plot of Synthesized Ehull vs. Normal Distribution")
plt.ylabel("Quantile", fontsize = 10)
plt.xlabel("$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize = 10)

#ax.text(
#    -0.195,    # x position, just left of the left spine
#    1.0,     # y position, just above the top spine
#    "B",      # the label
#    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
#    fontsize=11,              # size of the letter
#    fontweight="bold",        # make it bold
#    va="top",                 # vertical alignment
#    ha="left"                 # horizontal alignment
#)

plt.tight_layout()
plt.savefig("Figure3", dpi=1500, bbox_inches='tight')
plt.show()


