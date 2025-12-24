from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

df = pd.read_csv("data/external/elemental_reference_DFT_ESEN_UMA_12_25.csv")
OUTPUT_DIR = Path("figures/figures_s17")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

elements = df["element"].tolist()
x = np.arange(len(elements))
width = 0.38
FS = 26

fig, ax = plt.subplots(figsize=(25, 9))
ax.bar(
    x - width / 2, df["DFT_minus_UMA"].to_numpy(dtype=float), width, label="UMA-ODAC"
)
ax.bar(
    x + width / 2,
    df["DFT_minus_ESEN"].to_numpy(dtype=float),
    width,
    label="eSEN-ODAC25",
)

ax.axhline(0, linewidth=1.2)

ax.set_xticks(x, elements, rotation=90, fontsize=FS)
ax.set_ylabel(r"$\Delta{E}_{\text{total}}$ [DFT - MLIP] (eV/atom)", fontsize=FS)
ax.legend(fontsize=FS)

ax.yaxis.set_minor_locator(AutoMinorLocator(2))  # one minor between majors
ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
ax.tick_params(axis="y", which="major", labelsize=FS, length=8, width=1.2)
ax.tick_params(axis="y", which="minor", labelsize=FS - 2, length=4, width=1.0)

ax.set_axisbelow(True)
ax.grid(which="major", axis="y", linestyle="-", linewidth=1.0, alpha=0.65)
ax.grid(which="minor", axis="y", linestyle="-", linewidth=1.8, alpha=0.65)
ax.grid(which="major", axis="x", linestyle="-", linewidth=1.8, alpha=0.65)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "plot.png", dpi=600, bbox_inches="tight")
