from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

forme = [-0.492680, -0.4929992, -0.4929993, -0.49306676, -0.49306677]
force_tol = [0.02, 0.025, 0.03, 0.035, 0.04]

rel_forme = np.zeros(5)
for i in range(len(forme)):
    print(i)
    rel_forme[i] = forme[i] - forme[0]

# 5) Create scatter plot
fig, ax = plt.subplots(figsize=(3.8, 2.485))
plt.scatter(force_tol, rel_forme, s=15, alpha=1, edgecolors="none")
plt.xlabel(r"Force Convergence Criterion (eV/${\mathrm{\AA}}$)", fontsize=10)
plt.ylabel(
    "MOF-5 $ΔE_{\\mathrm{form}}$ w.r.t. Strictest \nConvergence (eV/atom)", fontsize=9
)


ax.tick_params(which="major", direction="in", length=10, width=1.25)
ax.tick_params(which="minor", direction="in", length=5, width=1.25)

for spine in ax.spines.values():
    spine.set_linewidth(1.25)

ax.tick_params(labelsize=8)

ax.minorticks_on()

# plt.legend(handles=[
#    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#6495ED', markersize=12, label='Synthesized'),
#    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6347', markersize=12, label='Not Synthesized')
# ],
# fontsize = 22,
# loc = 'best',
# frameon = False
# )
plt.xlim([0.017, 0.041])
plt.ylim([-0.0005, 0.0001])
# plt.ylim([-0.5, -0.470])

# ax.text(
#    -0.145,    # x position, just left of the left spine
#    1.0,     # y position, just above the top spine
#    "B",      # the label
#    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
#    fontsize=11,              # size of the letter
#    fontweight="bold",        # make it bold
#    va="top",                 # vertical alignment
#    ha="left"                 # horizontal alignment
# )

# plt.xlim([0, 40])
# plt.ylim([0, 0.82])

ax.text(
    -0.45,  # x position, just left of the left spine
    1.05,  # y position, just above the top spine
    "B",  # the label
    transform=ax.transAxes,  # interpret x,y in axis fraction (0 to 1)
    fontsize=11,  # size of the letter
    fontweight="bold",  # make it bold
    va="top",  # vertical alignment
    ha="left",  # horizontal alignment
)


plt.tight_layout()
plt.savefig("FigureS4b", dpi=1500, bbox_inches="tight")
#:wplt.show()
# s#ns.set_theme(style="ticks")
