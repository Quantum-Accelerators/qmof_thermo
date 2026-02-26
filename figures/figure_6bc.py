import json
import ast
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from pymatgen.core.periodic_table import Element
from matplotlib.lines import Line2D
import gzip

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Load data
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# ==================== DATA PREPARATION FOR PLOT B ====================
# HSAB Classification Dictionary
hsab_simple = {
    # Hard acids
    'Li_1': 'Hard', 'Na_1': 'Hard', 'K_1': 'Hard', 'Be_2': 'Hard', 'Mg_2': 'Hard', 'Ca_2': 'Hard', 'Sr_2': 'Hard',
    'Al_3': 'Hard', 'Ga_3': 'Hard', 'In_3': 'Hard', 'Sc_3': 'Hard', 'La_3': 'Hard', 'Ti_4': 'Hard',
    'Zr_4': 'Hard', 'Cr_3': 'Hard', 'Fe_3': 'Hard', 'Co_3': 'Hard', 'Th_4': 'Hard',
    'Pu_4': 'Hard', 'Yb_3': 'Hard', 'Mn_2': 'Hard', 'Gd_3': 'Hard', 'Mn_3': 'Hard', 'Sn_4': 'Hard', 'Lu_3': 'Hard',
    'As_3': 'Hard', 'Si_4': 'Hard', 'U_4': 'Hard', 'Ce_3': 'Hard', 'Hf_4': 'Hard', 'V_3': 'Hard', 'Mo_6': 'Hard', 'Np_6': 'Hard', 'Ce_4': 'Hard',
    # Soft acids
    'Cs_1': 'Soft', 'Tl_1': 'Soft', 'Tl_3': 'Soft', 'Pd_2': 'Soft', 'Pt_2': 'Soft', 'Cu_1': 'Soft', 'Ag_1': 'Soft',
    'Au_1': 'Soft', 'Cd_2': 'Soft', 'Hg_1': 'Soft', 'Hg_2': 'Soft', 'Pt_4': 'Soft', 'Te_4': 'Soft',
    # Borderline acids
    'Pb_2': 'Borderline', 'Sb_3': 'Borderline', 'Bi_3': 'Borderline', 'Fe_2': 'Borderline', 'Co_2': 'Borderline',
    'Ni_2': 'Borderline', 'Cu_2': 'Borderline', 'Zn_2': 'Borderline', 'Sn_2': 'Borderline', 'Ru_2': 'Borderline',
    'Rh_3': 'Borderline', 'Ir_3': 'Borderline', 'Os_2': 'Borderline', 'Rh_2':'Borderline', 'Mo_2': 'Borderline',
    'Cr_2': 'Borderline',
    # MOF Paper assignments
    'Rb_1': 'Hard', 'U_6': 'Hard', 'Y_3': 'Hard', 'Sm_3': 'Hard', 'Eu_3': 'Hard',
    'Ba_2': 'Hard', 'Tb_3': 'Hard', 'Pr_3': 'Hard', 'Nd_3': 'Hard',
    'Er_3': 'Hard', 'Tm_3': 'Hard', 'Yb': 'Hard', 'Dy_3': 'Hard', 'Ho_3': 'Hard'
}


hsab_colors = {
    'Hard': '#1f77b4',
    'Soft': '#d62728',
    'Borderline': '#9467bd',
}

def get_hsab_classification(element, oxidation_state):
    key = f"{element}_{oxidation_state}"
    return hsab_simple.get(key, 'Unknown')

def analyze_oxidation_states(oxi_data):
    if oxi_data is None:
        return False, None, None, "no_oximachine_data"
    metal_symbols = oxi_data.get('metal_symbols', [])
    predictions = oxi_data.get('prediction', [])
    max_probas = oxi_data.get('max_probas', [])
    if not metal_symbols or not predictions or not max_probas:
        return False, None, None, "incomplete_oximachine_data"
    unique_metals = set(metal_symbols)
    if len(unique_metals) != 1:
        return False, None, None, "multiple_metal_types"
    metal_type = list(unique_metals)[0]
    unique_oxidations = set(predictions)
    if len(unique_oxidations) != 1:
        return False, metal_type, None, "disagreeing_oxidation_states"
    oxidation_state = list(unique_oxidations)[0]
    if max(max_probas) < 0.85:
        return False, metal_type, oxidation_state, "low_confidence"
    return True, metal_type, oxidation_state, "valid"

def categorize_chemsys(chemsys):
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

# Process data for Plot B
valid_entries = []
unknown_hsab = []

# Track filtering reasons
total_mofs = len(data)
synthesized_count = sum(1 for d in data.values() if d.get('synthesizable', False))

filtering_stats = {
    'no_oximachine_data': 0,
    'incomplete_oximachine_data': 0,
    'multiple_metal_types': 0,
    'disagreeing_oxidation_states': 0,
    'low_confidence': 0,
    'not_N_or_O_only': 0,
    'unknown_hsab': 0,
    'valid': 0
}

for qmof_id, qmof_data in data.items():
    if not qmof_data.get('synthesizable', False):
        continue
    
    is_valid, metal_type, oxidation_state, reason = analyze_oxidation_states(qmof_data.get('oxidation_states'))
    
    if is_valid:
        chemsys = qmof_data.get('chemsys', '')
        no_category = categorize_chemsys(chemsys)
        if no_category in ["N_only", "O_only"]:
            hsab_class = get_hsab_classification(metal_type, oxidation_state)
            
            # Track unknown classifications
            if hsab_class == 'Unknown':
                filtering_stats['unknown_hsab'] += 1
                unknown_hsab.append({
                    'qmof_id': qmof_id,
                    'metal_oxi': f"{metal_type}_{oxidation_state}",
                    'no_category': no_category
                })
                continue
            
            filtering_stats['valid'] += 1
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
        else:
            filtering_stats['not_N_or_O_only'] += 1
    else:
        filtering_stats[reason] += 1

# Print comprehensive summary
print("\n" + "="*60)
print("PLOT B DATA FILTERING SUMMARY")
print("="*60)
print(f"Total MOFs in database: {total_mofs}")
print(f"Synthesized MOFs: {synthesized_count} ({synthesized_count/total_mofs*100:.1f}%)")
print(f"\nFiltering breakdown (synthesized only):")
print(f"  ✓ Valid for plotting: {filtering_stats['valid']} ({filtering_stats['valid']/synthesized_count*100:.1f}%)")
print(f"\nExcluded:")
print(f"  No oximachine data: {filtering_stats['no_oximachine_data']} ({filtering_stats['no_oximachine_data']/synthesized_count*100:.1f}%)")
print(f"  Incomplete oximachine data: {filtering_stats['incomplete_oximachine_data']} ({filtering_stats['incomplete_oximachine_data']/synthesized_count*100:.1f}%)")
print(f"  Multiple metal types: {filtering_stats['multiple_metal_types']} ({filtering_stats['multiple_metal_types']/synthesized_count*100:.1f}%)")
print(f"  Disagreeing oxidation states: {filtering_stats['disagreeing_oxidation_states']} ({filtering_stats['disagreeing_oxidation_states']/synthesized_count*100:.1f}%)")
print(f"  Low confidence (<85%): {filtering_stats['low_confidence']} ({filtering_stats['low_confidence']/synthesized_count*100:.1f}%)")
print(f"  Not N-only or O-only: {filtering_stats['not_N_or_O_only']} ({filtering_stats['not_N_or_O_only']/synthesized_count*100:.1f}%)")
print(f"  Unknown HSAB classification: {filtering_stats['unknown_hsab']} ({filtering_stats['unknown_hsab']/synthesized_count*100:.1f}%)")

total_excluded = synthesized_count - filtering_stats['valid']
print(f"\nTotal excluded: {total_excluded} ({total_excluded/synthesized_count*100:.1f}%)")

# Report unknown HSAB classifications
if unknown_hsab:
    df_unknown = pd.DataFrame(unknown_hsab)
    unknown_counts = df_unknown['metal_oxi'].value_counts()
    
    print(f"\n=== MISSING HSAB CLASSIFICATIONS ===")
    print(f"Metal-Oxidation combinations not in hsab_simple dictionary:")
    for combo, count in unknown_counts.items():
        print(f"  {combo}: {count} MOFs")

df_b = pd.DataFrame(valid_entries)
print("="*60 + "\n")
# ==================== DATA PREPARATION FOR PLOT C ====================
df_c = (
    pd.DataFrame.from_dict(data, orient="index")
      .reset_index()
      .rename(columns={"index": "qmof_id"})
)
df_c = df_c.loc[:, ~df_c.columns.duplicated()]
df_c = df_c[["qmof_id", "ehull", "linker_types", "synthesizable", "chemsys"]]
df_c = df_c[df_c["synthesizable"] == True].copy()

single_linkers = [
    "carboxylate_no_nitrogen",
    "imidazolate_no_oxygen",
    "pyrazine_no_oxygen",
   # "azide_no_oxygen",
   # "nitrile_no_oxygen",
    "triazole_no_oxygen",
    "tetrazine_no_oxygen",
    "pyrazole_no_oxygen",
    "tetrazolate_no_oxygen"
]

mixed_linker_combinations = []

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

df_c['linker_types'] = df_c['linker_types'].apply(ensure_list)

def categorize_linkers(linker_list):
    if not isinstance(linker_list, list) or len(linker_list) == 0:
        return None
    if (len(linker_list) == 2 and
        'tetrazine_no_oxygen' in linker_list and
        'triazole_no_oxygen' in linker_list):
        return 'triazole_no_oxygen'
    if len(linker_list) == 1 and linker_list[0] in single_linkers:
        return linker_list[0]
    if len(linker_list) > 1:
        sorted_linkers = tuple(sorted(linker_list))
        if sorted_linkers in mixed_linker_combinations:
            return " + ".join(sorted_linkers)
    return None

df_c['linker_category'] = df_c['linker_types'].apply(categorize_linkers)
df_c_filtered = df_c[df_c['linker_category'].notna()].copy()
df_c_filtered['linker_key'] = df_c_filtered['linker_category']

counts = df_c_filtered['linker_key'].value_counts()
median_ehull_c = df_c_filtered.groupby("linker_key")["ehull"].median().sort_values()
sorted_categories = median_ehull_c.index.tolist()

def simplify_linker_name(linker_name):
    simplified = linker_name.replace('_no_oxygen', '').replace('_no_nitrogen', '')
    simplified = simplified.replace('triazole', 'triazolate')
    simplified = simplified.replace('pyrazole', 'pyrazolate')
    if simplified == 'nitrile':
        simplified = 'nitrile + other'
    elif simplified == 'azide':
        simplified = 'azide + other'
    return simplified.replace('_', ' ').title()

# ==================== CREATE COMBINED FIGURE ====================
# Width = 3.4 + 3.4 = 6.8, Height = 2.75
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.85))

# ============== PLOT B: HSAB Ionic Radius vs Ehull ==============
# Prepare data for grouped violin plot
df_b['group'] = (
    df_b['no_category'].str.replace('_', '-') + '\n' + 
    df_b['hsab_class'].str.replace('_', '-')
)

# Define order
group_order = []
for element in ['O-only', 'N-only']:
    for hsab in ['Hard', 'Borderline', 'Soft']:
        group_name = f'{element}\n{hsab}'
        if group_name in df_b['group'].to_numpy():
            group_order.append(group_name)

print(f"Group order: {group_order}")

# Set categorical ordering
df_b['group'] = pd.Categorical(
    df_b['group'], 
    categories=group_order, 
    ordered=True
)

# Get counts for each group
group_counts = df_b['group'].value_counts()

# Get median ehull for each group
group_medians = df_b.groupby('group')['ehull'].median()
print(group_medians)



# Create red gradient (light to dark) for O-only
red_gradient = sns.color_palette("Reds", n_colors=5)[1:4]
red_gradient = red_gradient[::-1]
# Create blue gradient (light to dark) for N-only  
blue_gradient = sns.color_palette("Blues", n_colors=5)[1:4]
blue_gradient = blue_gradient[::-1]
# Combine into palette dictionary
palette = dict(zip(group_order, list(red_gradient) + list(blue_gradient)))

# Create violin plot
sns.violinplot(
    data=df_b,
    x='group',
    y='ehull',
    order=group_order,
    ax=ax1,
    inner='quartile',  # Shows quartiles
    color='lightgray',
    linewidth=1.25,
    palette = palette,
    cut=0
)

# Styling
ax1.tick_params(which='major', direction='in', length=10, width=1.25)
ax1.tick_params(which='minor', direction='in', length=5, width=1.25)
ax1.tick_params(labelsize=8)
ax1.tick_params(axis='y', pad=10)

for spine in ax1.spines.values():
    spine.set_linewidth(1.25)

ax1.set_ylabel('$ΔE_{\mathrm{hull}}$ (eV/atom)', fontsize=10)
ax1.set_xlabel('')
#ax1.set_ylim(0, 0.8)
ax1.set_xlim(-0.7, len(group_order) - 0.3)


# Set x-axis labels with counts
group_names = [
    'O-only\nHard', 
    'O-only\nBL', 
    'O-only\nSoft',
    'N-only\nHard', 
    'N-only\nBL', 
    'N-only\nSoft'
]

name_mapping = {
    'O-only\nHard': 'O-only\nHard',
    'O-only\nBL': 'O-only\nBorderline',
    'O-only\nSoft': 'O-only\nSoft',
    'N-only\nHard': 'N-only\nHard',
    'N-only\nBL': 'N-only\nBorderline',
    'N-only\nSoft': 'N-only\nSoft'
}

labels_with_counts = [
    f"{name}\n({group_counts.get(name_mapping[name], 0)})" 
    for name in group_names
]

ax1.set_xticklabels(labels_with_counts, fontsize=8, ha='center')
# Shift Borderline labels down
#labels = ax1.get_xticklabels()
#for i, label in enumerate(labels):
#    if "Borderline" in label.get_text():
#        label.set_y(-0.02)

ax1.set_ylim(-0.07,0.6)

ax1.text(-0.195, 1, "B", transform=ax1.transAxes, fontsize=11,
         fontweight="bold", va="top", ha="left")
# ============== PLOT C: Linker Type Violin Plots ==============
if len(sorted_categories) > 0:
    colors_c = sns.color_palette("Set1", len(sorted_categories))
    
    sns.violinplot(
        data=df_c_filtered,
        x="linker_key",
        y="ehull",
        order=sorted_categories,
        inner="quart",
        ax=ax2,
        width=0.9,
        cut=0,
        palette=colors_c,
        linewidth=0.9
    )
    
    labels_with_info = []
    for cat in sorted_categories:
        simplified_name = simplify_linker_name(cat)
        count = int(counts.get(cat, 0))
        label = f"{simplified_name}\n({count})"
        labels_with_info.append(label)
    
    ax2.set_xticklabels(labels_with_info, rotation=13.5, fontsize=8, ha='center')
    
    ax2.tick_params(which='major', direction='in', length=7, width=1.25)
    ax2.tick_params(which='minor', direction='in', length=5, width=1.25)
    ax2.tick_params(axis='y', which='both', length=0)
    ax2.tick_params(axis='y', labelsize=8)
    
    for spine in ax2.spines.values():
        spine.set_linewidth(1.25)
    
    ax2.set_ylabel("$ΔE_{\mathrm{hull}}$ (eV/atom)", fontsize=10)
    ax2.set_xlabel("")
    ax2.set_ylim([-0.03, 0.5])
    
    for y in ax2.get_yticks():
        ax2.axhline(y=y, linestyle="-", linewidth=1.05, color="grey", alpha=0.8, zorder=0)
    
    ax2.text(-0.16, 1.00, "C", transform=ax2.transAxes, fontsize=11,
             fontweight="bold", va="top", ha="left")

plt.tight_layout()
plt.savefig("Figure6_bc.png", dpi=1500, bbox_inches='tight')
plt.show()
