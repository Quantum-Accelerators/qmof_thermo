from __future__ import annotations

import gzip
import json

import matplotlib.pyplot as plt
import numpy as np

# Load both JSON files
with gzip.open("All_qmof_results_within_0.005_eVperatom.json.gz", "rt") as f:
    data_all = json.load(f)

with gzip.open("All_qmof_results_within_0.0_eVperatom.json.gz", "rt") as f:
    data_filtered = json.load(f)

# Calculate differences for each qmof_id
differences = []
qmof_ids = []
for qmof_id in data_all:
    if qmof_id in data_filtered:
        ehull_all = data_all[qmof_id].get("ehull")
        ehull_filtered = data_filtered[qmof_id].get("ehull")

        if ehull_all is not None and ehull_filtered is not None:
            diff = ehull_all - ehull_filtered
            differences.append(diff)
            qmof_ids.append(qmof_id)

print(f"Total entries with both ehull values: {len(differences)}")
print(f"Mean difference: {np.mean(differences):.6f} eV/atom")
print(f"Median difference: {np.median(differences):.6f} eV/atom")
print(f"Min difference: {np.min(differences):.6f} eV/atom")
print(f"Max difference: {np.max(differences):.6f} eV/atom")

filter_val = 0.002

# Split data into two groups
differences_below = [d for d in differences if d != 0 and d <= filter_val]
differences_above = [d for d in differences if d > filter_val]
n_no_change = sum(1 for d in differences if d == 0)

print(f"\nBelow/at filter ({filter_val}): {len(differences_below)}")
print(
    f"Above filter: {len(differences_above)} ({len(differences_above) / len(differences) * 100:.1f}%)"
)
print(f"No change: {n_no_change} ({n_no_change / len(differences) * 100:.1f}%)")

# Print outlier values
if len(differences_above) > 0:
    print(f"\nOutlier values (difference > {filter_val} eV/atom):")
    outliers = [
        (qmof_ids[i], differences[i])
        for i in range(len(differences))
        if differences[i] > filter_val
    ]
    outliers.sort(key=lambda x: x[1], reverse=True)
    for qmof_id, diff in outliers[:10]:  # Show top 10
        print(f"  {qmof_id}: {diff:.6f} eV/atom")
    if len(outliers) > 10:
        print(f"  ... and {len(outliers) - 10} more")

# Create side-by-side histogram
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Use same normalization base for both plots (all non-zero differences)
total_nonzero = len(differences_below) + len(differences_above)

# Left panel: data below/at filter
if len(differences_below) > 0:
    weights1 = np.ones_like(differences_below) / total_nonzero
    n1, bins1, patches1 = ax1.hist(
        differences_below,
        bins=12,
        weights=weights1,
        color="#6495ED",
        edgecolor="black",
        alpha=0.7,
    )
    ax1.set_xlabel(
        "Nonzero ΔΔ$E_{\\mathrm{hull}}$ (eV/atom)\nMP$_{\\mathrm{within}}$ $_{0.000}$ $_{\\mathrm{eV/atom}}$ → MP$_{\\mathrm{within}}$ $_{0.005}$ $_{\\mathrm{eV/atom}}$",
        fontsize=14,
    )
    ax1.set_ylabel("Normalized Frequency", fontsize=14)
    ax1.set_title(
        r"Main Distribution of nonzero ΔΔ$E_{\mathrm{hull}}$"
         f" (≤ {filter_val} eV/atom)\n$n$ = {len(differences_below)}",
        fontsize=14,
    )

    changed_total = len(differences_below) + len(differences_above)
    textstr = f"Nonzero ΔΔ$E_{{\\mathrm{{hull}}}}$: {changed_total}\nZero ΔΔ$E_{{\\mathrm{{hull}}}}$: {n_no_change}"
    ax1.text(
        0.95,
        0.95,
        textstr,
        transform=ax1.transAxes,
        fontsize=14,
        verticalalignment="top",
        horizontalalignment="right",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.8,
            "edgecolor": "black",
        },
    )
#    ax1.locator_params(axis='x', nbins=5)
# Right panel: data above filter
if len(differences_above) > 0:
    weights2 = np.ones_like(differences_above) / total_nonzero
    n2, bins2, patches2 = ax2.hist(
        differences_above,
        bins=12,
        weights=weights2,
        color="#FF6347",
        edgecolor="black",
        alpha=0.7,
    )
    ax2.set_xlabel(
        "Nonzero ΔΔ$E_{\\mathrm{hull}}$ (eV/atom)\nMP$_{\\mathrm{within}}$ $_{0.000}$ $_{\\mathrm{eV/atom}}$ → MP$_{\\mathrm{within}}$ $_{0.005}$ $_{\\mathrm{eV/atom}}$",
        fontsize=14,
    )
    ax2.set_ylabel("Normalized Frequency", fontsize=14)
    ax2.set_title(
        r"Outliers of nonzero ΔΔ$E_{\mathrm{hull}}$"
         f"(> {filter_val} eV/atom)\n$n$ = {len(differences_above)}",
        fontsize=14,
    )

# Styling for both axes
for ax in [ax1, ax2]:
    ax.minorticks_on()
    ax.tick_params(which="major", direction="in", length=8, width=1.5)
    ax.tick_params(which="minor", direction="in", length=4, width=1.2)
    ax.tick_params(labelsize=12)
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.grid(False)

ax1.text(
    -1.29,  # x position, just left of the left spine
    1.1,  # y position, just above the top spine
    "A",  # the label
    transform=ax.transAxes,  # interpret x,y in axis fraction (0 to 1)
    fontsize=24,  # size of the letter
    fontweight="bold",  # make it bold
    va="top",  # vertical alignment
    ha="left",  # horizontal alignment
)

ax2.text(
    -0.14,  # x position, just left of the left spine
    1.1,  # y position, just above the top spine
    "B",  # the label
    transform=ax.transAxes,  # interpret x,y in axis fraction (0 to 1)
    fontsize=24,  # size of the letter
    fontweight="bold",  # make it bold
    va="top",  # vertical alignment
    ha="left",  # horizontal alignment
)


plt.tight_layout()
plt.savefig("FigureS7ab.png", dpi=300, bbox_inches="tight")
plt.show()

# Print statistics
positive_changes = sum(1 for d in differences if d > 0)
negative_changes = sum(1 for d in differences if d < 0)
no_change = sum(1 for d in differences if d == 0)
changed = positive_changes + negative_changes

print("\nDistribution of changes:")
print(f"  Changed: {changed} ({changed / len(differences) * 100:.1f}%)")
print(f"  No change: {no_change} ({no_change / len(differences) * 100:.1f}%)")
print("\nOf those that changed:")
print(
    f"  Positive (ehull increased): {positive_changes} ({positive_changes / len(differences) * 100:.1f}%)"
)
print(
    f"  Negative (ehull decreased): {negative_changes} ({negative_changes / len(differences) * 100:.1f}%)"
)

threshold = 0.01
above_threshold = sum(1 for d in differences if d > threshold)
print(
    f"\nEntries with difference > {threshold} eV/atom: {above_threshold} ({above_threshold / len(differences) * 100:.1f}%)"
)
