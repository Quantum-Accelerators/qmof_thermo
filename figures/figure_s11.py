from __future__ import annotations

import gzip
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import auc, roc_curve

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

df = df[["ehull", "synthesizable"]].copy()

# Assuming df has 'ehull' and 'synthesizable' columns
y_true = df["synthesizable"].astype(int)  # True=1, False=0
y_scores = -df["ehull"]  # Negative because lower is better

fpr, tpr, thresholds = roc_curve(y_true, y_scores)
actual_thresholds = -thresholds
roc_auc = auc(fpr, tpr)

optimal_idx = np.argmax(tpr - fpr)
print(np.argmax(tpr - fpr))
print(tpr[np.argmax(tpr - fpr)])
print(fpr[np.argmax(tpr - fpr)])

optimal_ehull = -thresholds[optimal_idx]

plt.figure(figsize=(10, 8))
plt.plot(fpr, tpr, "b-", lw=2, label=f"ROC (AUC = {roc_auc:.3f})")
plt.plot([0, 1], [0, 1], "k--", lw=1, label="Random")

# Select specific thresholds to highlight
threshold_values = [
    optimal_ehull,
    0.267,
    0.306,
    0.353,
    0.476,
]  # Ehull values you care about
colors = ["red", "orange", "green", "purple", "blue"]

for thresh, color in zip(threshold_values, colors, strict=False):
    # Find closest actual threshold
    idx = np.argmin(np.abs(actual_thresholds - thresh))
    actual_val = actual_thresholds[idx]
    actual_val = round(actual_val, 3)
    # Plot point
    plt.plot(
        fpr[idx],
        tpr[idx],
        "o",
        color=color,
        markersize=11,
        label=r"Δ$E_{\mathrm{hull}}$ = " + str(actual_val) + " eV/atom",
    )

    # Add text label next to point
    plt.text(
        fpr[idx],
        tpr[idx] - 0.1,
        f"{actual_val:.3f}\n(TPR={tpr[idx]:.2f})",
        fontsize=12,
        color=color,
        fontweight="bold",
        ha="center",
        va="bottom",
    )

plt.xlabel("False Positive Rate for Hypothetical MOFs", fontsize=14)
plt.ylabel("True Positive Rate for Synthesized MOFs", fontsize=14)
plt.title(r"ROC Curve with Δ$E_{\mathrm{hull}}$ Thresholds", fontsize=16)
plt.legend(loc="lower right", fontsize=11, frameon=True)
plt.grid(alpha=0.3)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.tight_layout()
plt.savefig("FigureS11", dpi=1000)
plt.show()
