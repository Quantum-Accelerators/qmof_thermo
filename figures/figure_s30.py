import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gzip

# Load synthesizable MOFs only
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

df = (
    pd.DataFrame.from_dict(data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)

# Filter for synthesized only
df_syn = df[df['synthesizable'] == True][['formation_energy', 'ehull']].copy()
df_syn['synthesizable'] = True

# Load new hypothetical MOFs
with gzip.open("ARC_MOFs.json.gz", "rt") as f:
    data2 = json.load(f)

df2 = (
    pd.DataFrame.from_dict(data2, orient="index")
      .reset_index()
      .rename(columns={"index": "id"})
)
df_hypo = df2[['formation_energy', 'ehull']].copy()
df_hypo['synthesizable'] = False

# Combine datasets
df_combined = pd.concat([df_syn, df_hypo], ignore_index=True)

print(f"Synthesized: {len(df_syn)}")
print(f"Hypothetical (New): {len(df_hypo)}")

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.6, 2.55))

# ============== PLOT A (Formation Energy) ==============
df_plot1 = df_combined[['formation_energy', 'synthesizable']].copy()
df_plot1['group'] = 'All MOFs'

sns.violinplot(
    data=df_plot1,
    x='group',
    y='formation_energy',
    hue='synthesizable',
    hue_order=[True, False],
    split=True,
    inner='quart',
    palette={True: '#6495ED', False: '#90EE90'},
    ax=ax1
)

ax1.set_xticklabels([])
ax1.tick_params(which='major', direction='in', length=10, width=1.25)
ax1.tick_params(which='minor', direction='in', length=5, width=1.25)
ax1.tick_params(axis='y', labelsize=9)
ax1.tick_params(axis='x', which='both', length=0)
ax1.minorticks_on()
ax1.set_ylim([-1.5, 0.65])
ax1.set_xlabel('', fontsize=1)
ax1.set_ylabel('Δ$E_{\mathrm{form}}$ (eV/atom)', fontsize=10)

for spine in ax1.spines.values():
    spine.set_linewidth(1.25)

ax1.text(-0.265, 1.00, 'A', transform=ax1.transAxes, fontsize=11,
         fontweight='bold', va='top', ha='left')
ax1.get_legend().remove()


# ============== PLOT B (E_hull) ==============
df_plot2 = df_combined[['ehull', 'synthesizable']].copy()
df_plot2['group'] = 'All MOFs'

print(df_plot2['synthesizable'].value_counts())

sns.violinplot(
    data=df_plot2,
    x='group',
    y='ehull',
    hue='synthesizable',
    hue_order=[True, False],
    split=True,
    inner='quart',
    palette={True: '#6495ED', False: '#90EE90'},
    ax=ax2
)

ax2.set_xticklabels([])
ax2.tick_params(which='major', direction='in', length=10, width=1.25)
ax2.tick_params(which='minor', direction='in', length=5, width=1.25)
ax2.tick_params(axis='y', labelsize=9)
ax2.tick_params(axis='x', which='both', length=0)
ax2.minorticks_on()
ax2.set_ylim([0, 0.8])
ax2.set_xlabel('', fontsize=22)
ax2.set_ylabel('$ΔE_{\mathrm{hull}}$ (eV/atom)', fontsize=10)

for spine in ax2.spines.values():
    spine.set_linewidth(1.25)

# Update legend labels
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, ['Syn QMOF', 'ARC-MOF'], 
           title='', fontsize=8,
           frameon=False, loc='upper right',
           bbox_to_anchor=([1.02, 1]))

ax2.text(-0.22, 1.00, 'B', transform=ax2.transAxes, fontsize=11,
         fontweight='bold', va='top', ha='left')

plt.tight_layout()
plt.savefig('Figure1.png', dpi=1500, bbox_inches='tight')
plt.show()
