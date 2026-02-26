import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core.periodic_table import Element
from adjustText import adjust_text
import gzip
from matplotlib.lines import Line2D


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# HSAB Classification Dictionary
hsab_simple = {
    # Hard acids
    'Li_1': 'Hard', 'Na_1': 'Hard', 'K_1': 'Hard', 'Be_2': 'Hard', 'Mg_2': 'Hard', 'Ca_2': 'Hard', 'Sr_2': 'Hard',
    'Al_3': 'Hard', 'Ga_3': 'Hard', 'In_3': 'Hard', 'Sc_3': 'Hard', 'La_3': 'Hard', 'Ti_4': 'Hard',
    'Zr_4': 'Hard', 'Cr_3': 'Hard', 'Fe_3': 'Hard', 'Co_3': 'Hard', 'Co_3': 'Hard', 'Th_4': 'Hard',
    'Pu_4': 'Hard', 'Yb_3': 'Hard', 'Mn_2': 'Hard', 'Gd_3': 'Hard', 'Mn_3': 'Hard', 'Sn_4': 'Hard', 'Lu_3': 'Hard',
    'As_3': 'Hard', 'Si_4': 'Hard', 'Zr_4': 'Hard', 'U_4': 'Hard', 'Ce_3': 'Hard', 'Hf_4': 'Hard', 
    
    # Soft acids
    'Cs_1': 'Soft', 'Tl_1': 'Soft', 'Tl_3': 'Soft', 'Pd_2': 'Soft', 'Pt_2': 'Soft', 'Cu_1': 'Soft', 'Ag_1': 'Soft',
    'Au_1': 'Soft', 'Cd_2': 'Soft', 'Hg_1': 'Soft', 'Hg_2': 'Soft', 'Pt_4': 'Soft', 'Te_4': 'Soft', 
    
    # Borderline acids
    'Pb_2': 'Borderline', 'Sb_3': 'Borderline', 'Bi_3': 'Borderline', 'Fe_2': 'Borderline', 'Co_2': 'Borderline',
    'Ni_2': 'Borderline', 'Cu_2': 'Borderline', 'Zn_2': 'Borderline', 'Sn_2': 'Borderline', 'Ru_2': 'Borderline',
    'Rh_3': 'Borderline', 'Ir_3': 'Borderline', 'Os_2': 'Borderline',

    # MOF Paper assignments, Cation exhange in MOFs: the HSAB principle appraisal
    'Rb_1': 'Hard', 'U_6': 'Hard', #from UO2 (2+)
    'Y_3': 'Hard', 'Sm_3': 'Hard', 'Eu_3': 'Hard',

    # "Recent Advances in the Study of Trivalent Lanthanides..
    'Ba_2': 'Hard', 'Tb_3': 'Hard', 'Pr_3': 'Hard', 'Nd_3': 'Hard', 
    'Er_3': 'Hard', 'Tm_3': 'Hard', 'Yb': 'Hard', 'Dy_3': 'Hard',  'Ho_3': 'Hard'

}

# HSAB Color mapping
hsab_colors = {
    'Hard': '#1f77b4',      # Blue
    'Soft': '#d62728',      # Red  
    'Borderline': '#9467bd', # Purple
 #   'Unknown': '#7f7f7f'    # Gray
}

def get_hsab_classification(element, oxidation_state):
    """Get HSAB classification using the simple dictionary lookup"""
    key = f"{element}_{oxidation_state}"
    return hsab_simple.get(key, 'Unknown')

# Load your results JSON with oximachine data
with gzip.open("All_qmof_results.json.gz", "r") as f:
    data = json.load(f)

def analyze_oxidation_states(oxi_data):
    """
    Analyze oxidation state data and return filtering results
    Returns: (is_valid, metal_type, oxidation_state, reason)
    """
    if oxi_data is None:
        return False, None, None, "no_oximachine_data"

    metal_symbols = oxi_data.get('metal_symbols', [])
    predictions = oxi_data.get('prediction', [])
    max_probas = oxi_data.get('max_probas', [])

    if not metal_symbols or not predictions or not max_probas:
        return False, None, None, "incomplete_oximachine_data"

    # Check if only one type of metal
    unique_metals = set(metal_symbols)
    if len(unique_metals) != 1:
        return False, None, None, "multiple_metal_types"

    metal_type = list(unique_metals)[0]

    # Check if all oxidation states agree
    unique_oxidations = set(predictions)
    if len(unique_oxidations) != 1:
        return False, metal_type, None, "disagreeing_oxidation_states"

    oxidation_state = list(unique_oxidations)[0]

    # Check confidence level (at least one prediction >= 0.95)
    if max(max_probas) < 0.85:
        return False, metal_type, oxidation_state, "low_confidence"

    return True, metal_type, oxidation_state, "valid"

def categorize_chemsys(chemsys):
    """Categorize MOF based on presence of N and O in chemsys"""
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

# Initialize data structures
valid_entries = []
filtering_stats = {
    "no_oximachine_data": 0,
    "incomplete_oximachine_data": 0,
    "multiple_metal_types": 0,
    "disagreeing_oxidation_states": 0,
    "low_confidence": 0,
    "valid": 0
}

# Store examples for each exception case
examples = {
    "no_oximachine_data": [],
    "incomplete_oximachine_data": [],
    "multiple_metal_types": [],
    "disagreeing_oxidation_states": [],
    "low_confidence": []
}

# Additional filtering stats for N/O categories
no_filtering_stats = {
    "N_and_O": 0,
    "neither": 0,
    "N_only": 0,
    "O_only": 0
}

# Filter MOFs based on criteria
for qmof_id, qmof_data in data.items():
    if not qmof_data.get('synthesizable', False):
        continue

    is_valid, metal_type, oxidation_state, reason = analyze_oxidation_states(qmof_data.get('oxidation_states'))
    filtering_stats[reason] += 1

    # Store examples (up to 3 per category)
    if reason != "valid" and len(examples[reason]) < 3:
        examples[reason].append(qmof_id)

    if is_valid:
        chemsys = qmof_data.get('chemsys', '')
        no_category = categorize_chemsys(chemsys)
        no_filtering_stats[no_category] += 1
        
        # Only keep N-only or O-only MOFs
        if no_category in ["N_only", "O_only"]:
            # Get HSAB classification
            hsab_class = get_hsab_classification(metal_type, oxidation_state)
            
            valid_entries.append({
                'qmof_id': qmof_id,
                'ehull': qmof_data['ehull'],
                'metal_type': metal_type,
                'oxidation_state': oxidation_state,
                'metal_oxi_combo': f"{metal_type}_{oxidation_state}",
                'chemsys': chemsys,
                'no_category': no_category,
                'metal_oxi_no_combo': f"{metal_type}_{oxidation_state}_{no_category}",
                'hsab_class': hsab_class
            })

print("Filtering Statistics:")
print(f"Total synthesizable MOFs: {sum(filtering_stats.values())}")
for reason, count in filtering_stats.items():
    print(f"  {reason}: {count}")
    if reason != "valid" and examples[reason]:
        print(f"    Examples: {', '.join(examples[reason])}")

print("\nN/O Category Statistics:")
for category, count in no_filtering_stats.items():
    print(f"  {category}: {count}")

print(f"\nExcluded N+O and neither: {no_filtering_stats['N_and_O'] + no_filtering_stats['neither']}")
print(f"Kept N-only and O-only: {no_filtering_stats['N_only'] + no_filtering_stats['O_only']}")

print()

# Create DataFrame from valid entries
df = pd.DataFrame(valid_entries)

if df.empty:
    print("No valid entries found after filtering!")
    exit()

# Group by metal-oxidation-N/O combination and calculate median ehull
median_ehull_by_group = df.groupby('metal_oxi_no_combo').agg({
    'ehull': 'median',
    'metal_type': 'first',
    'oxidation_state': 'first',
    'no_category': 'first',
    'hsab_class': 'first',
    'qmof_id': 'count'  # Count entries per group
}).reset_index()

# Rename count column for clarity
median_ehull_by_group = median_ehull_by_group.rename(columns={'qmof_id': 'count'})

# Filter out groups with less than 8 entries
median_ehull_by_group = median_ehull_by_group[median_ehull_by_group['count'] >= 8]

print(f"Groups after filtering (≥8 entries): {len(median_ehull_by_group)}")

def get_ionic_radius(metal_type, oxidation_state):
    """Get ionic radius using pymatgen"""

    element = Element(metal_type)
    ionic_radius = element.ionic_radii.get(oxidation_state, None)
    return ionic_radius
    
# Add ionic radius to the dataframe
median_ehull_by_group['ionic_radius'] = median_ehull_by_group.apply(
    lambda row: get_ionic_radius(row['metal_type'], row['oxidation_state']), axis=1
)

# Remove any rows where ionic radius couldn't be found
median_ehull_by_group = median_ehull_by_group.dropna(subset=['ionic_radius'])

print(f"Final dataset: {len(median_ehull_by_group)} metal-oxidation-N/O combinations")

# Print HSAB distribution
print("\nHSAB Classification Distribution:")
hsab_counts = median_ehull_by_group['hsab_class'].value_counts()
for hsab_class, count in hsab_counts.items():
    print(f"  {hsab_class}: {count}")

fig, ax = plt.subplots(figsize=(15, 13))

# Separate data by N/O category and plot with different markers, colored by HSAB
o_only_data = median_ehull_by_group[median_ehull_by_group['no_category'] == 'O_only']
n_only_data = median_ehull_by_group[median_ehull_by_group['no_category'] == 'N_only']

# Plot O-only data with circles, colored by HSAB
for hsab_class, color in hsab_colors.items():
    subset = o_only_data[o_only_data['hsab_class'] == hsab_class]
    if not subset.empty:
        ax.scatter(subset["ionic_radius"], subset["ehull"],
                  s=140, c=color, alpha=1,
                  marker='o',
                  label=f"{hsab_class} (O-only)",
                  edgecolors='black', linewidth=1)

# Plot N-only data with triangles
for hsab_class, color in hsab_colors.items():
    subset = n_only_data[n_only_data['hsab_class'] == hsab_class]
    if not subset.empty:
        ax.scatter(subset["ionic_radius"], subset["ehull"],
                  s=140, c=color, alpha=1,
                  marker='^',
                  label=f"{hsab_class} (N-only)",
                  edgecolors='black', linewidth=1)
# Add element labels with oxidation states
#texts = []
#for idx, row in median_ehull_by_group.iterrows():
#    label = f"{row['metal_type']}(+{row['oxidation_state']})"
#    text = ax.annotate(label,
#                      (row['ionic_radius'], row['ehull']),
#                      fontsize=12, ha='left', va='bottom',
#                      color='black')
#    texts.append(text)

# Adjust text positions to avoid overlap
#adjust_text(texts, expand_points=(2, 1), force_points=2, iter_lim=100, max_move=3)

# Styling
ax.tick_params(which='major', direction='in', length=26, width=1.8)
ax.tick_params(which='minor', direction='in', length=14, width=1.6)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(labelsize=20)
ax.minorticks_on()
ax.tick_params(axis='x', which='major', pad=12)
ax.tick_params(axis='y', which='major', pad=12)

# Add labels
plt.xlabel('Ionic Radius (Å)', fontsize=22)
plt.ylabel('Median Δ$E_{\mathrm{hull}}$ (eV/atom)', fontsize=22)


# Legend for HSAB classification (colors)
hsab_legend_elements = [Line2D([0], [0], marker='o', color='w', 
                              markerfacecolor=hsab_colors[hsab_class], markersize=10,
                              markeredgecolor='black', markeredgewidth=1,
                              label=hsab_class) for hsab_class in sorted(hsab_colors.keys())]

# Legend for ligand types (shapes)
ligand_legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=10,
           markeredgecolor='black', markeredgewidth=1, label='O-only ligands'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', markersize=10,
           markeredgecolor='black', markeredgewidth=1, label='N-only ligands')
]

# Add both legends
hsab_legend = ax.legend(handles=hsab_legend_elements, title='Metal HSAB Classification', 
                       loc='upper right', frameon=False, fontsize=20, title_fontsize=20)
ax.add_artist(hsab_legend)  # Keep the first legend when adding the second

ligand_legend = ax.legend(handles=ligand_legend_elements, title='Ligand Type', bbox_to_anchor=(0.95, 0.775), 
                         loc='center right', frameon=False, fontsize=20, title_fontsize=20)

#plt.xlim([0.5,1.85])

# Add panel label
#ax.text(-0.15, 1.045, "B", transform=ax.transAxes, fontsize=24, 
#        fontweight="bold", va="top", ha="left")

plt.tight_layout()
plt.savefig("FigureS19", dpi=1000, bbox_inches='tight')
plt.show()

# Print enhanced summary data with HSAB classification
print("\nSummary of metal-oxidation-ligand combinations with HSAB classification:")
print("=" * 130)
for idx, row in median_ehull_by_group.iterrows():
    ligand_type = "O-only" if row['no_category'] == 'O_only' else "N-only"
    marker = "○" if row['no_category'] == 'O_only' else "△"
    print(f"{row['metal_type']}(+{row['oxidation_state']}) {ligand_type} {marker}: ehull={row['ehull']:.4f}, ionic_radius={row['ionic_radius']:.3f}, count={row['count']}, HSAB={row['hsab_class']}")

print(f"\nTotal data points plotted: {len(median_ehull_by_group)}")
print(f"O-only groups: {len(o_only_data)}")
print(f"N-only groups: {len(n_only_data)}")

# HSAB analysis summary
print(f"\nHSAB Analysis Summary:")
for hsab_class in ['Hard', 'Borderline', 'Soft', 'Unknown']:
    subset = median_ehull_by_group[median_ehull_by_group['hsab_class'] == hsab_class]
    if not subset.empty:
        avg_ehull = subset['ehull'].mean()
        print(f"  {hsab_class} acids: {len(subset)} groups, avg ehull = {avg_ehull:.4f} eV/atom")
