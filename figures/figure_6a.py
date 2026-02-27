from __future__ import annotations

import gzip
import json

import matplotlib.pyplot as plt
import pandas as pd
from periodic_trends import plotter

pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", None)
pd.set_option("display.max_colwidth", None)

# 1) Load your results JSON
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

metal_list = [
    "Ag",
    "Al",
    "Au",
    "Ba",
    "Be",
    "Bi",
    "Ca",
    "Cd",
    "Ce",
    "Co",
    "Cr",
    "Cs",
    "Cu",
    "Dy",
    "Er",
    "Eu",
    "Fe",
    "Ga",
    "Gd",
    "Hf",
    "Hg",
    "Ho",
    "In",
    "Ir",
    "K",
    "La",
    "Li",
    "Lu",
    "Mg",
    "Mn",
    "Mo",
    "Na",
    "Nb",
    "Nd",
    "Ni",
    "Np",
    "Pb",
    "Pd",
    "Pr",
    "Pt",
    "Pu",
    "Rb",
    "Re",
    "Rh",
    "Ru",
    "Sc",
    "Sm",
    "Sn",
    "Sr",
    "Tb",
    "Tc",
    "Th",
    "Ti",
    "Tl",
    "Tm",
    "U",
    "V",
    "W",
    "Y",
    "Yb",
    "Zn",
    "Zr",
    "B",
    "As",
    "Sb",
    "At",
    "Te",
    "Si",
    "Ge",
]


def has_single_metal(frac_dict):
    """Returns True if MOF has exactly one type of metal"""
    metals_present = [metal for metal in metal_list if metal in frac_dict]
    return len(metals_present) == 1


def get_single_metal(frac_dict):
    """Returns the single metal if MOF has exactly one type, None otherwise"""
    metals_present = [metal for metal in metal_list if metal in frac_dict]
    return metals_present[0] if len(metals_present) == 1 else None


# Filter data to only include single-metal MOFs
single_metal_data = {}
multi_metal_count = 0
no_metal_count = 0

for qmof_id, properties in data.items():
    if not properties.get("synthesizable", False):
        continue

    frac_comp = properties.get("frac_composition", {})
    metals_present = [metal for metal in metal_list if metal in frac_comp]

    if len(metals_present) == 1:
        single_metal_data[qmof_id] = properties
    elif len(metals_present) > 1:
        multi_metal_count += 1
    else:
        no_metal_count += 1

print("Filtering results:")
print(f"  Single-metal MOFs: {len(single_metal_data)}")
print(f"  Multi-metal MOFs (excluded): {multi_metal_count}")
print(f"  No-metal MOFs (excluded): {no_metal_count}")

# 2) Build a tidy DataFrame from filtered data
df = (
    pd.DataFrame.from_dict(single_metal_data, orient="index")
    .reset_index()
    .rename(columns={"index": "qmof_id"})
)

df = df[["ehull", "synthesizable", "chemsys", "frac_composition"]]

# Create expanded dataframe - each MOF will have exactly one metal
expanded_rows = []
for _idx, row in df.iterrows():
    metal = get_single_metal(row["frac_composition"])
    if metal:  # Should always be true due to filtering
        new_row = row.copy()
        new_row["group"] = metal
        expanded_rows.append(new_row)

df = pd.DataFrame(expanded_rows).reset_index(drop=True)

print("\nMOF distribution by metal (single-metal only):")
metal_counts = df["group"].value_counts()
print(metal_counts)

# Filter out metals with fewer than 10 entries
metals_with_enough_data = metal_counts[metal_counts >= 8].index
df_filtered = df[df["group"].isin(metals_with_enough_data)]

excluded_metals = metal_counts[metal_counts < 8]
print("\nExcluded metals with <10 entries:")
for metal, count in excluded_metals.items():
    print(f"  {metal}: {count} entries")

print(
    f"\nFinal dataset: {len(df_filtered)} MOFs with {len(metals_with_enough_data)} metals (≥10 entries each)"
)

median_ehull_by_metal = df_filtered.groupby("group")["ehull"].median().reset_index()
median_ehull_by_metal.columns = ["Element", "Median_Ehull"]

# Sort by median ehull
median_ehull_by_metal = median_ehull_by_metal.sort_values("Median_Ehull")

print("\nMedian Ehull by Metal/Metalloid (single-metal MOFs, ≥10 entries):")
print(median_ehull_by_metal)

# Save to CSV
median_ehull_by_metal.to_csv(
    "single_metal_median_ehull_filtered.csv", index=False, header=False
)
print("\nSaved results to 'single_metal_median_ehull_filtered.csv'")

# Optional: Also create a version with additional statistics
detailed_stats = (
    df_filtered.groupby("group")["ehull"]
    .agg(["count", "median", "mean", "std", "min", "max"])
    .reset_index()
)
detailed_stats.columns = [
    "Element",
    "Count",
    "Median_Ehull",
    "Mean_Ehull",
    "Std_Ehull",
    "Min_Ehull",
    "Max_Ehull",
]
detailed_stats = detailed_stats.sort_values("Median_Ehull")
detailed_stats.to_csv("single_metal_detailed_ehull_stats_filtered.csv", index=False)
print(
    "Also saved detailed statistics to 'single_metal_detailed_ehull_stats_filtered.csv'"
)

fig, ax = plt.subplots(figsize=(6.5, 3))
plotter("single_metal_median_ehull_filtered.csv")
plt.savefig("Figure6a")

# plt.show()
