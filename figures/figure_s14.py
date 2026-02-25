import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gzip
# Load data
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# Build DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)
#df = df[["ehull", "pore_diameters", "synthesizable", "chemsys"]]
#df["pore_diameters"] = df["pore_diameters"].str[0]

print("=== Data Coverage ===")
print(f"Total entries: {len(df)}")
print(f"Entries with pore_diameters: {df['pore_diameters'].notna().sum()}")
print(f"Entries WITHOUT pore_diameters (skipped): {df['pore_diameters'].isna().sum()}")
print(f"Percentage with data: {df['pore_diameters'].notna().sum()/len(df)*100:.1f}%")

missing_pores = df[df['pore_diameters'].isna()]['qmof_id'].values.tolist()
if missing_pores:
    print(f"\n=== Entries WITHOUT pore_diameters ({len(missing_pores)}) ===")
    for qmof_id in missing_pores[:20]:  # Show first 20
        print(f"  {qmof_id}")
    if len(missing_pores) > 20:
        print(f"  ... and {len(missing_pores) - 20} more")
else:
    print("\nAll entries have pore diameter data!")

df = df[["formation_energy", "ehull", "pore_diameters", "synthesizable", "chemsys"]]
df["pore_diameters"] = df["pore_diameters"].str[0]


print("All data counts:")
print(df["synthesizable"].value_counts())

# Create 1x2 subplot grid
# Width = 3.4 + 3.4 = 6.8, Height = max(2.5, 2.485) = 2.5
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.5))

# ============== PLOT A: All MOFs (Synthesized + Unsynthesized) ==============
color_map = {True: "#6495ED", False: "#FF6347"}
colors_all = df["synthesizable"].map(color_map)

ax1.scatter(df["pore_diameters"], df["ehull"], c=colors_all, s=4, alpha=0.4, edgecolors='none')
ax1.set_xlabel("Largest Cavity Diameter (Å)", fontsize=10)
ax1.set_ylabel("$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=10)
ax1.tick_params(which='major', direction='in', length=10, width=1.25)
ax1.tick_params(which='minor', direction='in', length=5, width=1.25)
ax1.tick_params(labelsize=8)
ax1.minorticks_on()

for spine in ax1.spines.values():
    spine.set_linewidth(1.25)


ax1.legend(handles=[
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#6495ED', markersize=5, label='Synthesized'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6347', markersize=5, label='Hypothetical')
], fontsize=9, loc='upper right', frameon=False, handletextpad=0.1)


ax1.set_xlim([0, 46])
ax1.set_ylim([0, 0.82])

ax1.text(-0.145, 1.0, "A", transform=ax1.transAxes, fontsize=11,
         fontweight="bold", va="top", ha="left")

# ============== PLOT B: Only Synthesized MOFs ==============
#df_syn = df[df["synthesizable"] == True]
#print("\nSynthesized data counts:")
#print(df_syn["synthesizable"].value_counts())

#colors_syn = df_syn["synthesizable"].map(color_map)

ax2.scatter(df["pore_diameters"], df["formation_energy"], c=colors_all, s=4, alpha=0.4, edgecolors='none')
ax2.set_xlabel("Largest Cavity Diameter (Å)", fontsize=10)
ax2.set_ylabel("$ΔE_{\mathrm{form}}$ (eV/atom)", fontsize=10)
ax2.tick_params(which='major', direction='in', length=10, width=1.25)
ax2.tick_params(which='minor', direction='in', length=5, width=1.25)
ax2.tick_params(labelsize=8)
ax2.minorticks_on()

for spine in ax2.spines.values():
    spine.set_linewidth(1.25)

#ax2.set_xlim([0, 38])
#ax2.set_ylim([0, 0.82])

ax2.text(-0.145, 1.0, "B", transform=ax2.transAxes, fontsize=11,
         fontweight="bold", va="top", ha="left")

plt.tight_layout()
plt.savefig("FigureS14.png", dpi=1500, bbox_inches='tight')
plt.show()
