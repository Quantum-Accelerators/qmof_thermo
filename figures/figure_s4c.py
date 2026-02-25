import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

kpt_Zn = [0.10003897,	0.15003897,	0.20003897,	0.25003897,	0.300038970]
kpt_O2 = [0.10486093,	0.15486093,	0.20486093,	0.25486093,	0.304860928]
kpt_C = [0.10005411,	0.15005411,	0.20005411,	0.25005411,	0.300054106 ]
kpt_H2 = [0.35000000,	0.40000000,	0.45000000,	0.50000000,	0.549999996]

total_nrg_Zn =[-1.523146565,	-1.522932245,	-1.52353486,	-1.520353485,	-1.528828045]
total_nrg_O2=[-5.000332055,	-5.000329354,	-5.000326885,	-5.000335748,	-5.000335766]
total_nrg_C = [-9.36641969,	-9.366546305,	-9.36645497,	-9.36674973,	-9.36681791]
total_nrg_H2=[-3.390882745,	-3.390882745,	-3.3908986,	-3.3908986,	-3.3908986]


#kspacing = ["BTK - 0.1", "BTK - 0.05", "BTK", "BTK + 0.05", "BTK + 0.1"]
forme=[-0.493036054,	-0.492995557,	-0.492988451,	-0.493092905,	-0.492422437]

rel_forme=np.zeros(5)
for i in range(0,len(forme)):
    print(i)
    rel_forme[i] = forme[i] - forme[0]


kspacing = ["0.1", "0.05", "0", "-0.05", "-0.1"]

kspacing = kspacing[::-1]
#form_e_MOF_5 = form_e_MOF_5[::-1]


# 5) Create scatter plot
fig, ax = plt.subplots(figsize=(4.6, 2.485))
plt.scatter(kspacing, rel_forme, s=15, alpha=1, edgecolors='none')
plt.xlabel(r"Shift in KSPACING Value (${\mathrm{\AA^{-1}}}$)", fontsize = 10 )
plt.ylabel("MOF-5 $ΔE_{\mathrm{form}}$ w.r.t.Strictest\n Convergence (eV/atom)", fontsize = 9)


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
# Categories are at positions 0, 1, 2, 3, 4
plt.xlim([-0.5, 4.5])  # Add padding on both sides
# Or show only first 3
plt.ylim([-0.0002, 0.0007])
#plt.xlim([350, 700])
#plt.ylim([-0.5, -0.470])

ax.text(
    -0.32,    # x position, just left of the left spine
    1.05,     # y position, just above the top spine
    "C",      # the label
    transform=ax.transAxes,   # interpret x,y in axis fraction (0 to 1)
    fontsize=11,              # size of the letter
    fontweight="bold",        # make it bold
    va="top",                 # vertical alignment
    ha="left"                 # horizontal alignment
)


plt.tight_layout()
plt.savefig("FigureS4C", dpi=1500, bbox_inches='tight')
#plt.show()


