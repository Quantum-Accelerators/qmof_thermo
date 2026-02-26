import json
import ast
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import gzip

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Load the combined JSON with both linker and oxidation state data
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# Build DataFrame
df = pd.DataFrame.from_dict(data, orient="index").reset_index().rename(columns={"index": "qmof_id"})

# Filter for synthesizable MOFs only
df = df[df["synthesizable"] == True].copy()

# Helper function to ensure linker_types is a list
def ensure_list(x):
    if x is None:
        return []
    if isinstance(x, (list, tuple)):
        return list(x)
    if isinstance(x, str):
        try:
            parsed = ast.literal_eval(x)
            if isinstance(parsed, (list, tuple)):
                return list(parsed)
            return [parsed]
        except Exception:
            return [x]
    return [x]

df['linker_types'] = df['linker_types'].apply(ensure_list)

def categorize_chemsys(chemsys):
    """Categorize MOF based on presence of N and O in chemsys"""
    if pd.isna(chemsys):
        return "unknown"
    elements = chemsys.split('-')
    has_nitrogen = 'N' in elements
    has_oxygen = 'O' in elements

    if has_oxygen and not has_nitrogen:
        return "O_only"
    elif has_nitrogen and not has_oxygen:
        return "N_only"
    elif has_nitrogen and has_oxygen:
        return "N_and_O"
    else:
        return "neither"

# Add N/O categorization
df['no_category'] = df['chemsys'].apply(categorize_chemsys)

# Filter to keep only N_only or O_only MOFs
df = df[df['no_category'].isin(['N_only', 'O_only'])].copy()

# Define azolate linkers
azolate_linkers = {
    "tetrazolate_no_oxygen",
    "triazole_no_oxygen", 
    "imidazolate_no_oxygen",
    "pyrazole_no_oxygen"
}

# Define target metal-oxidation combinations
target_metals = {
    'Zn_2': 'Zn²⁺',
#    'Co_2': 'Co²⁺',
    'Cd_2': 'Cd²⁺',
    'Ag_1': 'Ag⁺',
    'Cu_2': 'Cu²⁺',
#    'Cu_1': 'Cu⁺'
 #   'Mn_2': 'Mn²⁺',
 #   'Mg_2': 'Mg²⁺',
 #   'Hg_2': 'Hg²⁺',


}

def categorize_linkers(linker_list):
    """
    Categorize linkers as:
    - 'Carboxylate': only carboxylate_no_nitrogen
    - 'Azolate': only azolate linkers (can be multiple types)
    - None: anything else
    """
    if not isinstance(linker_list, list) or len(linker_list) == 0:
        return None
    
    linker_set = set(linker_list)
    
    # Check for carboxylate-only
    if linker_set == {"carboxylate_no_nitrogen"}:
        return "Carboxylate"
    
    # Check for azolate-only (can be multiple azolate types)
    if linker_set.issubset(azolate_linkers) and len(linker_set) > 0:
        return "Azolate"
    
    return None

# Apply linker categorization
df['linker_category'] = df['linker_types'].apply(categorize_linkers)

# Filter to keep only categorized linkers
df_filtered = df[df['linker_category'].notna()].copy()

# Filter for MOFs with monometallic cations
df_filtered = df_filtered[df_filtered['monometallic_cation'].notna()].copy()

# Create metal-oxidation key
def create_metal_key(mono_cation):
    if mono_cation is None:
        return None
    return f"{mono_cation['element']}_{mono_cation['oxidation_state']}"

df_filtered['metal_key'] = df_filtered['monometallic_cation'].apply(create_metal_key)

# Filter for target metals only
df_filtered = df_filtered[df_filtered['metal_key'].isin(target_metals.keys())].copy()

# Create display labels with counts
metal_counts = {}
for metal_key, metal_label in target_metals.items():
    metal_data = df_filtered[df_filtered['metal_key'] == metal_key]
    carboxylate_count = len(metal_data[metal_data['linker_category'] == 'Carboxylate'])
    azolate_count = len(metal_data[metal_data['linker_category'] == 'Azolate'])
    label_with_counts = f"{metal_label}\n(C: {carboxylate_count}, A: {azolate_count})"
    metal_counts[metal_key] = label_with_counts

df_filtered['metal_label'] = df_filtered['metal_key'].map(metal_counts)

print(f"Data after filtering: {len(df_filtered)} MOFs")
print("\nN/O Category distribution:")
print(df_filtered['no_category'].value_counts())
print("\nMetal distribution:")
print(df_filtered['metal_key'].value_counts())
print("\nLinker category distribution:")
print(df_filtered['linker_category'].value_counts())
print("\nMetal-Linker combinations:")
print(pd.crosstab(df_filtered['metal_key'], df_filtered['linker_category']))

# Create the split violin plot
fig, ax = plt.subplots(figsize=(12, 8))

# Define colors for each linker type
colors = ['#ff7f0e', '#2ca02c']  # Orange for Carboxylate, Green for Azolate

# Create split violin plot
sns.violinplot(
    data=df_filtered,
    x="metal_label", 
    y="ehull",
    hue="linker_category",
    split=True,
    inner="quart",
    palette=colors,
    ax=ax,
    cut=0
)

# Styling
ax.tick_params(which='major', direction='in', length=26, width=1.8)
ax.tick_params(which='minor', direction='in', length=14, width=1.6)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(labelsize=20)

ax.minorticks_on()
ax.tick_params(axis='x', which='major', pad=12)
ax.tick_params(axis='y', which='major', pad=12)

# Labels and title
plt.ylabel("Δ$E_{\mathrm{hull}}$ (eV/atom)", fontsize=22)
plt.xlabel("Metal Cation", fontsize=22)

# Customize legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, title='Linker Type', fontsize=16, title_fontsize=18, 
          frameon=False, loc='upper right')

# Add grid for better readability
#ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.8)

ax.tick_params(axis='x', which='both', length=0)

# Set y-axis limits
plt.ylim([0.0, 0.5])

plt.tight_layout()
plt.savefig("FigureS21.png", dpi=1200, bbox_inches='tight')
plt.show()

# Print detailed statistics
print("\nDetailed Statistics by Metal and Linker Type (N/O filtered):")
print("=" * 80)

for metal_key, metal_label in target_metals.items():
    print(f"\n{metal_label}:")
    metal_data = df_filtered[df_filtered['metal_key'] == metal_key]
    for linker_cat in ['Carboxylate', 'Azolate']:
        subset = metal_data[metal_data['linker_category'] == linker_cat]
        if len(subset) > 0:
            median_ehull = subset['ehull'].median()
            mean_ehull = subset['ehull'].mean()
            std_ehull = subset['ehull'].std()
            count = len(subset)
            # Show N/O breakdown
            no_breakdown = subset['no_category'].value_counts()
            print(f"  {linker_cat}: {count} MOFs")
            print(f"    Median: {median_ehull:.4f} eV/atom")
            print(f"    Mean ± Std: {mean_ehull:.4f} ± {std_ehull:.4f} eV/atom")
            print(f"    N/O breakdown: {dict(no_breakdown)}")
        else:
            print(f"  {linker_cat}: 0 MOFs")

# Show some example linker combinations for azolates
print("\nExample azolate combinations:")
azolate_data = df_filtered[df_filtered['linker_category'] == 'Azolate']
linker_combinations = azolate_data['linker_types'].apply(lambda x: ', '.join(sorted(x))).value_counts()
print(linker_combinations.head(10))

# Statistical comparison
print("\nStatistical Summary:")
summary_stats = df_filtered.groupby(['metal_key', 'linker_category'])['ehull'].agg([
    'count', 'median', 'mean', 'std'
]).round(4)
print(summary_stats)
