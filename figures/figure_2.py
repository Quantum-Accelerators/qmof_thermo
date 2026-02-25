import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import numpy as np
import gzip

shiftx = -0.02
shifty = 0.03
# --- Load all data files ---
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    qmof_data = json.load(f)


with gzip.open("MP_data_with_MOF_chemsys.json.gz", "rt") as f:
    mp_data = json.load(f)

# Fix 1: Use mp_data, not data
mp_df = (
    pd.DataFrame.from_dict(mp_data, orient="index")  # Use mp_data here
      .reset_index()
      .rename(columns={"index": "material_id"})
)

# Fix 2: Correct filtering syntax
# Fix 3: .values (no parentheses)
mp_syn_ehull = mp_df[mp_df["theoretical"] == True]["energy_above_hull"].values

mp_nonsyn_ehull = mp_df[mp_df["theoretical"] == False]["energy_above_hull"].values

mp_syn_eform = mp_df[mp_df["theoretical"] == True]["formation_energy_per_atom"].values

mp_nonsyn_eform = mp_df[mp_df["theoretical"] == False]["formation_energy_per_atom"].values

mp_ehull = np.concatenate([mp_syn_ehull, mp_nonsyn_ehull])
mp_eform = np.concatenate([mp_syn_eform, mp_nonsyn_eform])

# Create 2x2 subplot grid
# Width = 3.25 + 3.25 = 6.5, Height = 2.5 + 2.25 = 4.75
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(6.5, 4.75))

print("qmof data is " + str(len(qmof_data)))
print("MP syn is " + str(len(mp_syn_eform)))
print("MP syn is " + str(len(mp_syn_ehull)))


print("MP nonsyn is "+ str(len(mp_nonsyn_eform)))
print("MP nonsyn is "+ str(len(mp_nonsyn_ehull)))



# ============== PLOT A: QMOF Hexbin ==============
df_qmof = (
    pd.DataFrame.from_dict(qmof_data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)[["ehull", "formation_energy"]]

df_qmof = df_qmof[
    (df_qmof.formation_energy >= -1.3) & (df_qmof.formation_energy <= 0.4) &
    (df_qmof.ehull >= 0) & (df_qmof.ehull <= 1.5)
]

hb1 = ax1.hexbin(df_qmof["formation_energy"], df_qmof["ehull"],
                 gridsize=80, cmap='viridis', mincnt=1)

cb1 = plt.colorbar(hb1, ax=ax1, pad=0.02)
cb1.ax.tick_params(labelsize=9)
ax1.tick_params(axis='both', labelsize=9)
ax1.minorticks_on()
ax1.tick_params(which='major', direction='in', length=10, width=1.25)
ax1.tick_params(which='minor', direction='in', length=5, width=1.25)

for spine in ax1.spines.values():
    spine.set_linewidth(1.25)

ax1.set_xlabel("$ΔE_{\mathrm{form}}$ (eV/atom)", fontsize=10)
ax1.set_ylabel("$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=10)
ax1.set_title("MOFs in QMOF Database", fontsize=10)

ax1.text(-0.19+shiftx, 1.07+shifty, "A", transform=ax1.transAxes, fontsize=11,
         fontweight="bold", va="top", ha="left")

# ============== PLOT B: MP Hexbin ==============
df_mp = pd.DataFrame({"ehull": mp_ehull, "eform": mp_eform})
df_mp = df_mp[
    (df_mp.eform >= -6) & (df_mp.eform <= 2) &
    (df_mp.ehull >= 0) & (df_mp.ehull <= 0.81)
]

hb2 = ax2.hexbin(df_mp["eform"], df_mp["ehull"],
                 gridsize=30, cmap='viridis', mincnt=1)

cb2 = plt.colorbar(hb2, ax=ax2, pad=0.02)
cb2.ax.tick_params(labelsize=9)
ax2.tick_params(axis='both', labelsize=9)
ax2.minorticks_on()
ax2.tick_params(which='major', direction='in', length=10, width=1.25)
ax2.tick_params(which='minor', direction='in', length=5, width=1.25)

for spine in ax2.spines.values():
    spine.set_linewidth(1.25)

ax2.set_xlabel("$ΔE_{\mathrm{form}}$ (eV/atom)", fontsize=10)
ax2.set_ylabel("$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=10)
ax2.set_title("MP within MOF Chemical Systems", fontsize=10)

ax2.text(-0.2+shiftx, 1.07+shifty, "B", transform=ax2.transAxes, fontsize=11,
         fontweight="bold", va="top", ha="left")

# ============== PLOT C: Formation Energy Histogram ==============
mp_df_syn_c = pd.DataFrame({"formation_energy": mp_syn_eform})
mp_df_syn_c["Dataset"] = "MP within QMOF Chemsys (Synthesized)"

mp_df_nonsyn_c = pd.DataFrame({"formation_energy": mp_nonsyn_eform})
mp_df_nonsyn_c["Dataset"] = "MP within QMOF Chemsys (Hypothetical)"

qmof_df_c = (
    pd.DataFrame.from_dict(qmof_data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)[["formation_energy", "synthesizable", "chemsys"]]

qmof_df_syn_c = qmof_df_c[qmof_df_c["synthesizable"] == True][["formation_energy"]].assign(
    Dataset="All MOFs in QMOF (Synthesized)")
qmof_df_nonsyn_c = qmof_df_c[qmof_df_c["synthesizable"] == False][["formation_energy"]].assign(
    Dataset="All MOFs in QMOF (Hypothetical)")

df_c = pd.concat([mp_df_syn_c, mp_df_nonsyn_c, qmof_df_syn_c, qmof_df_nonsyn_c], ignore_index=True)

x_min_c, x_max_c = -4, 1.8
bin_width_c = 0.2
bin_edges_c = np.arange(x_min_c, x_max_c + bin_width_c, bin_width_c)

sns.histplot(data=df_c[df_c.Dataset == "MP within QMOF Chemsys (Synthesized)"],
             x="formation_energy", bins=bin_edges_c, stat="probability", element="step",
             color="#3CB371", fill=False, alpha=1, linewidth=1.25,
             label="MP within QMOF Chemsys (Synthesized)", ax=ax3)

sns.histplot(data=df_c[df_c.Dataset == "All MOFs in QMOF (Synthesized)"],
             x="formation_energy", bins=bin_edges_c, stat="probability", element="step",
             fill=False, color="#6495ED", alpha=1, linewidth=1.25,
             label="All MOFs in QMOF (Synthesized)", ax=ax3)

sns.histplot(data=df_c[df_c.Dataset == "All MOFs in QMOF (Hypothetical)"],
             x="formation_energy", bins=bin_edges_c, stat="probability", element="step",
             fill=False, color="#FF6347", alpha=1, linewidth=1.25,
             label="All MOFs in QMOF (Unsynthesized)", ax=ax3)

sns.histplot(data=df_c[df_c.Dataset == "MP within QMOF Chemsys (Hypothetical)"],
             x="formation_energy", bins=bin_edges_c, stat="probability", element="step",
             color="#FF8C00", fill=False, alpha=1, linewidth=1.25,
             label="MP within QMOF Chemsys (Unsynthesized)", ax=ax3)

ax3.tick_params(which='major', direction='in', length=12, width=1.25)
ax3.tick_params(which='minor', direction='in', length=6, width=1.25)
ax3.tick_params(labelsize=8)
ax3.minorticks_on()
ax3.tick_params(axis="x", pad=5)

for spine in ax3.spines.values():
    spine.set_linewidth(1.25)

ax3.set_xlabel("$ΔE_{\mathrm{form}}$ (eV/atom)", fontsize=10)
ax3.set_ylabel("Relative Frequency", fontsize=10)
ax3.set_xlim([x_min_c, x_max_c])

ax3.text(-0.15+shiftx, 1+shifty, "C", transform=ax3.transAxes, fontsize=11,
         fontweight="bold", va="top", ha="left")

# ============== PLOT D: Ehull Histogram ==============
mp_df_syn_d = pd.DataFrame({"ehull": mp_syn_ehull})
mp_df_syn_d["Dataset"] = "MP within QMOF Chemsys (Synthesized)"

mp_df_nonsyn_d = pd.DataFrame({"ehull": mp_nonsyn_ehull})
mp_df_nonsyn_d["Dataset"] = "MP within QMOF Chemsys (Hypothetical)"

print("this is the median for unsynthesized MP ehull")
print(mp_df_nonsyn_d["ehull"].median())

print("this is the median for synthesized MP ehull")
print(mp_df_syn_d["ehull"].median())

mp_df_syn_d = mp_df_syn_d[mp_df_syn_d["ehull"] < 0.51]
mp_df_nonsyn_d = mp_df_nonsyn_d[mp_df_nonsyn_d["ehull"] < 0.51]



qmof_df_d = (
    pd.DataFrame.from_dict(qmof_data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)[["ehull", "synthesizable", "chemsys"]]

qmof_df_syn_d = qmof_df_d[qmof_df_d["synthesizable"]][["ehull"]].assign(
    Dataset="All MOFs in QMOF (Synthesized)")
qmof_df_nonsyn_d = qmof_df_d[~qmof_df_d["synthesizable"]][["ehull"]].assign(
    Dataset="All MOFs in QMOF (Hypothetical)")

df_d = pd.concat([mp_df_syn_d, mp_df_nonsyn_d, qmof_df_syn_d, qmof_df_nonsyn_d], ignore_index=True)

x_min_d, x_max_d = 0, 0.5
bin_width_d = 0.02
bin_edges_d = np.arange(x_min_d, x_max_d + bin_width_d, bin_width_d)

sns.histplot(data=df_d[df_d.Dataset == "MP within QMOF Chemsys (Synthesized)"],
             x="ehull", bins=bin_edges_d, stat="probability", element="step",
             fill=False, color="#3CB371", linewidth=1.25, label="MP - Synthesized", ax=ax4)

sns.histplot(data=df_d[df_d.Dataset == "MP within QMOF Chemsys (Hypothetical)"],
             x="ehull", bins=bin_edges_d, stat="probability", element="step",
             fill=False, color="#FF8C00", linewidth=1.25, label="MP - Hypothetical", ax=ax4)

sns.histplot(data=df_d[df_d.Dataset == "All MOFs in QMOF (Synthesized)"],
             x="ehull", bins=bin_edges_d, stat="probability", element="step",
             fill=False, color="#6495ED", linewidth=1.25, label="QMOF - Synthesized", ax=ax4)

sns.histplot(data=df_d[df_d.Dataset == "All MOFs in QMOF (Hypothetical)"],
             x="ehull", bins=bin_edges_d, stat="probability", element="step",
             fill=False, color="#FF6347", linewidth=1.25, label="QMOF - Hypothetical", ax=ax4)

ax4.tick_params(which='major', direction='in', length=12, width=1.25)
ax4.tick_params(which='minor', direction='in', length=6, width=1.25)
ax4.tick_params(labelsize=8)
ax4.minorticks_on()
ax4.tick_params(axis="x", pad=5)

for spine in ax4.spines.values():
    spine.set_linewidth(1.25)

ax4.axvline(0, linestyle="--", color="gray")
ax4.set_xlabel("$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=10)
ax4.set_ylabel("Relative Frequency", fontsize=10)
ax4.set_xlim([x_min_d, x_max_d])
ax4.set_ylim([0, 0.56])

ax4.legend(title="Dataset", fontsize=8, title_fontsize=8, frameon=False, loc='upper right')

ax4.text(-0.15+shiftx, 1+shifty, "D", transform=ax4.transAxes, fontsize=11,
         fontweight="bold", va="top", ha="left")

plt.tight_layout()
plt.savefig("Figure2.png", dpi=1500, bbox_inches='tight')
plt.show()
