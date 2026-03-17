from __future__ import annotations

import ast
import gzip
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)

# HSAB classification dictionary
hsab_simple = {
    # Hard acids
    "Li_1": "Hard",
    "Na_1": "Hard",
    "K_1": "Hard",
    "Be_2": "Hard",
    "Mg_2": "Hard",
    "Ca_2": "Hard",
    "Sr_2": "Hard",
    "Al_3": "Hard",
    "Ga_3": "Hard",
    "In_3": "Hard",
    "Sc_3": "Hard",
    "La_3": "Hard",
    "Ti_4": "Hard",
    "Zr_4": "Hard",
    "Cr_3": "Hard",
    "Fe_3": "Hard",
    "Co_3": "Hard",
    "Th_4": "Hard",
    "Pu_4": "Hard",
    "Yb_3": "Hard",
    "Mn_2": "Hard",
    "Gd_3": "Hard",
    "Mn_3": "Hard",
    "Sn_4": "Hard",
    "Lu_3": "Hard",
    "As_3": "Hard",
    "Si_4": "Hard",
    "U_4": "Hard",
    "Ce_3": "Hard",
    "Hf_4": "Hard",
    # Soft acids
    "Cs_1": "Soft",
    "Tl_1": "Soft",
    "Tl_3": "Soft",
    "Pd_2": "Soft",
    "Pt_2": "Soft",
    "Cu_1": "Soft",
    "Ag_1": "Soft",
    "Au_1": "Soft",
    "Cd_2": "Soft",
    "Hg_1": "Soft",
    "Hg_2": "Soft",
    "Pt_4": "Soft",
    "Te_4": "Soft",
    # Borderline acids
    "Pb_2": "Borderline",
    "Sb_3": "Borderline",
    "Bi_3": "Borderline",
    "Fe_2": "Borderline",
    "Co_2": "Borderline",
    "Ni_2": "Borderline",
    "Cu_2": "Borderline",
    "Zn_2": "Borderline",
    "Sn_2": "Borderline",
    "Ru_2": "Borderline",
    "Rh_3": "Borderline",
    "Ir_3": "Borderline",
    "Os_2": "Borderline",
    # MOF Paper assignments
    "Rb_1": "Hard",
    "U_6": "Hard",
    "Y_3": "Hard",
    "Sm_3": "Hard",
    "Eu_3": "Hard",
    "Ba_2": "Hard",
    "Tb_3": "Hard",
    "Pr_3": "Hard",
    "Nd_3": "Hard",
    "Er_3": "Hard",
    "Tm_3": "Hard",
    "Yb": "Hard",
    "Dy_3": "Hard",
    "Ho_3": "Hard",
}

# HSAB Color mapping
hsab_colors = {
    "Hard": "#1f77b4",  # Blue
    "Soft": "#d62728",  # Red
    "Borderline": "#9467bd",  # Purple
    "Unknown": "#7f7f7f",  # Gray
}


def get_hsab_classification(element, oxidation_state):
    """Get HSAB classification using the simple dictionary lookup"""
    if pd.isna(element) or pd.isna(oxidation_state):
        return "Unknown"
    # Ensure oxidation_state is an integer
    try:
        oxidation_state = int(oxidation_state)
    except (ValueError, TypeError):
        return "Unknown"
    key = f"{element}_{oxidation_state}"
    return hsab_simple.get(key, "Unknown")


# 1) Load your results JSON
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# 2) Build a tidy DataFrame
df = (
    pd.DataFrame.from_dict(data, orient="index")
    .reset_index()
    .rename(columns={"index": "qmof_id"})
)

# Remove any duplicate columns
df = df.loc[:, ~df.columns.duplicated()]

# Keep the columns you need + monometallic_cation
df = df[
    [
        "qmof_id",
        "ehull",
        "linker_types",
        "synthesizable",
        "chemsys",
        "monometallic_cation",
    ]
]

df = df[df["synthesizable"]].copy()

print(df["synthesizable"].value_counts())
print("\nlinker_types")
print(df["linker_types"].value_counts())

# --- Define allowed single linkers and mixed combinations ---
single_linkers = [
    "carboxylate_no_nitrogen",
    "imidazolate_no_oxygen",
    "pyrazine_no_oxygen",
    "triazole_no_oxygen",
    "tetrazine_no_oxygen",
    "pyrazole_no_oxygen",
    "tetrazolate_no_oxygen",
]

mixed_linker_combinations = []


# Helper: ensure linker_types is a real Python list
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


df["linker_types"] = df["linker_types"].apply(ensure_list)


def categorize_linkers(linker_list):
    """Categorize a linker list"""
    if not isinstance(linker_list, list) or len(linker_list) == 0:
        return None

    # Special case: collapse tetrazine_no_oxygen + triazole_no_oxygen -> triazole_no_oxygen
    if (
        len(linker_list) == 2
        and "tetrazine_no_oxygen" in linker_list
        and "triazole_no_oxygen" in linker_list
    ):
        return "triazole_no_oxygen"

    # Check for single linkers
    if len(linker_list) == 1 and linker_list[0] in single_linkers:
        return linker_list[0]

    # Check for mixed combinations
    if len(linker_list) > 1:
        sorted_linkers = tuple(sorted(linker_list))
        if sorted_linkers in mixed_linker_combinations:
            return " + ".join(sorted_linkers)

    return None


# Apply categorization
df["linker_category"] = df["linker_types"].apply(categorize_linkers)

# Filter to keep only rows with recognized categories
df_filtered = df[df["linker_category"].notna()].copy()

print(f"Rows before filter: {len(df)}")
print(f"Rows after category filter: {len(df_filtered)}")


# Extract metal element and oxidation state from monometallic_cation
def extract_metal_info(cation_dict):
    """Extract the element and oxidation state from monometallic_cation dictionary"""
    if pd.isna(cation_dict) or not isinstance(cation_dict, dict):
        return None, None
    element = cation_dict.get("element", None)
    oxidation_state = cation_dict.get("oxidation_state", None)
    return element, oxidation_state


df_filtered[["metal_element", "oxidation_state"]] = df_filtered[
    "monometallic_cation"
].apply(lambda x: pd.Series(extract_metal_info(x)))

# Remove rows where metal_element is None
df_filtered = df_filtered[df_filtered["metal_element"].notna()].copy()

print(f"Rows after removing None metal elements: {len(df_filtered)}")

# Apply HSAB classification
df_filtered["hsab_class"] = df_filtered.apply(
    lambda row: get_hsab_classification(row["metal_element"], row["oxidation_state"]),
    axis=1,
)

# Print HSAB distribution
print("\nHSAB Classification Distribution:")
print(df_filtered["hsab_class"].value_counts())

# Debug: Print some examples of metal/oxidation state combinations
print("\n" + "=" * 60)
print("DEBUGGING: Sample metal/oxidation state combinations:")
print("=" * 60)
sample_df = df_filtered[
    ["qmof_id", "metal_element", "oxidation_state", "hsab_class"]
].head(20)
for _idx, row in sample_df.iterrows():
    key = f"{row['metal_element']}_{int(row['oxidation_state'])}"
    print(
        f"{row['qmof_id']}: {row['metal_element']} with ox state {row['oxidation_state']} -> key: {key} -> {row['hsab_class']}"
    )

# Print which keys are actually in the dictionary
print("\nSample keys in hsab_simple dictionary:")
print(list(hsab_simple.keys())[:20])


def simplify_linker_name(linker_name):
    """Simplify linker names by removing suffixes"""
    simplified = linker_name.replace("_no_oxygen", "").replace("_no_nitrogen", "")
    simplified = simplified.replace("triazole", "triazolate")
    simplified = simplified.replace("pyrazole", "pyrazolate")
    return simplified.replace("_", " ").title()


# Calculate proportions of each HSAB class for each linker category
hsab_counts = (
    df_filtered.groupby(["linker_category", "hsab_class"])
    .size()
    .reset_index(name="count")
)

# Calculate total MOFs per linker category
category_totals = (
    df_filtered.groupby("linker_category").size().reset_index(name="total")
)

# Merge to get proportions
hsab_counts = hsab_counts.merge(category_totals, on="linker_category")
hsab_counts["proportion"] = (hsab_counts["count"] / hsab_counts["total"]) * 100

# Sort categories by median ehull (to maintain consistent ordering)
median_ehull = df_filtered.groupby("linker_category")["ehull"].median().sort_values()
sorted_categories = median_ehull.index.tolist()

# Create pivot table for stacked bar plot
pivot_data = hsab_counts.pivot_table(
    index="linker_category", columns="hsab_class", values="proportion"
)
pivot_data = pivot_data.reindex(sorted_categories)  # Order by median ehull
pivot_data = pivot_data.fillna(0)  # Fill missing values with 0

# Get category counts for labels
counts = df_filtered["linker_category"].value_counts()

# Create the stacked bar plot
fig, ax = plt.subplots(figsize=(12, 8))

# Define order of HSAB classes for stacking (bottom to top)
hsab_order = ["Hard", "Borderline", "Soft", "Unknown"]
# Only include classes that exist in the data
hsab_order = [h for h in hsab_order if h in pivot_data.columns]

# Create the stacked bar plot
bottom = np.zeros(len(sorted_categories))
bars = []
for hsab_class in hsab_order:
    if hsab_class in pivot_data.columns:
        values = pivot_data[hsab_class].to_numpy()
        bar = ax.bar(
            range(len(sorted_categories)),
            values,
            bottom=bottom,
            label=hsab_class,
            color=hsab_colors[hsab_class],
            width=0.7,
        )
        bars.append(bar)
        bottom += values

# Customize the plot
labels_with_info = []
for cat in sorted_categories:
    simplified_name = simplify_linker_name(cat)
    count = int(counts.get(cat, 0))
    label = f"{simplified_name} ({count})"
    labels_with_info.append(label)

ax.set_xticks(range(len(sorted_categories)))
ax.set_xticklabels(labels_with_info, rotation=13.5, fontsize=20)

ax.tick_params(which="major", direction="in", length=26, width=1.8)
ax.tick_params(which="minor", direction="in", length=14, width=1.6)
for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(axis="y", which="both", length=0)
ax.tick_params(axis="y", labelsize=20)
plt.ylabel("Proportion of HSAB Class (%)", fontsize=22)
plt.xlabel("")
plt.ylim([0, 100])

# Add horizontal gridlines
for y in range(0, 101, 20):
    ax.axhline(y=y, linestyle="-", linewidth=1.3, color="grey", zorder=0, alpha=0.3)

# Add legend
ax.legend(
    title="HSAB Class",
    loc="upper right",
    bbox_to_anchor=(1.3, 1.05),
    fontsize=16,
    title_fontsize=16,
    frameon=False,
)

# Add label in top-left corner

ax.tick_params(axis="x", which="both", length=0)

plt.tight_layout()
plt.savefig("FigureS20.png", dpi=1500, bbox_inches="tight")
plt.show()

# Print summary statistics
print("\n" + "=" * 60)
print("SUMMARY STATISTICS:")
print("=" * 60)
print(f"Total synthesizable MOFs: {len(df)}")
print(f"MOFs in recognized categories: {len(df_filtered)}")

print("\nHSAB distribution by linker category:")
for category in sorted_categories:
    print(f"\n{simplify_linker_name(category)}:")
    cat_data = df_filtered[df_filtered["linker_category"] == category]
    hsab_dist = cat_data["hsab_class"].value_counts()
    total = len(cat_data)
    for hsab_class in ["Hard", "Borderline", "Soft", "Unknown"]:
        if hsab_class in hsab_dist.index:
            count = hsab_dist[hsab_class]
            percentage = (count / total) * 100
            print(f"  {hsab_class}: {count} ({percentage:.1f}%)")


print("METAL ION COUNTS BY LINKER CATEGORY:")
print("=" * 60)
for category in sorted_categories:
    print(f"\n{simplify_linker_name(category)} ({counts.get(category, 0)} total MOFs):")
    cat_data = df_filtered[df_filtered["linker_category"] == category]

    # Create metal ion identifier (element + oxidation state)
    cat_data_copy = cat_data.copy()
    cat_data_copy["metal_ion"] = cat_data_copy.apply(
        lambda row: f"{row['metal_element']}({int(row['oxidation_state'])}+)", axis=1
    )

    # Count each metal ion type
    metal_ion_counts = cat_data_copy["metal_ion"].value_counts()

    # Print top metal ions
    for metal_ion, count in metal_ion_counts.items():
        percentage = (count / len(cat_data)) * 100
        # Get HSAB class for this metal ion
        element, ox_str = metal_ion.split("(")
        ox_state = int(ox_str.replace("+)", ""))
        hsab = get_hsab_classification(element, ox_state)
        print(f"  {metal_ion}: {count} ({percentage:.1f}%) - {hsab}")
