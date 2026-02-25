import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

encut = [400, 520, 600, 680]
forme = [-0.478126657, -0.4930, -0.4928, -0.4929]

rel_forme=np.zeros(4)
for i in range(0,len(forme)):
    print(i)
    rel_forme[i] = forme[i] - forme[3]

print(rel_forme)
# 5) Create scatter plot
fig, ax = plt.subplots(figsize=(3.8, 2.485))
plt.scatter(encut, rel_forme, s=15, alpha=1, edgecolors='none')
plt.xlabel("ENCUT (eV)", fontsize = 10 )
plt.ylabel("MOF-5 $ΔE_{\mathrm{form}}$ w.r.t. \nStrictest Convergence (eV/atom)", fontsize = 9)


ax.tick_params(
    which='major', direction='in', length=10, width=1.25
)
ax.tick_params(
    which='minor', direction='in', length=5, width=1.25
)

for spine in ax.spines.values():
    spine.set_linewidth(1.25)

ax.tick_params(labelsize=8)

ax.minorticks_on()

#plt.legend(handles=[
#    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#6495ED', markersize=12, label='Synthesized'),
#    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6347', markersize=12, label='Not Synthesized')
#],
#fontsize = 22,
#loc = 'best',
#frameon = False
#)
plt.xlim([350, 700])
#plt.ylim([-0.5, -0.470])

ax.text(
    -0.37,    # x position, just left of the left spine
    1.05,     # y position, just above the top spine
    "A",      # the label
    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
    fontsize=11,              # size of the letter
    fontweight="bold",        # make it bold
    va="top",                 # vertical alignment
    ha="left"                 # horizontal alignment
)

#plt.xlim([0, 40])
plt.ylim([-0.002, 0.0155])

plt.tight_layout()
plt.savefig("FigureS4a", dpi=1500, bbox_inches='tight')
print("done")
#plt.show()
#s#ns.set_theme(style="ticks")


