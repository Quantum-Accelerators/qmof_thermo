import json
from collections import Counter
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
from ase import Atoms
import gzip

# Define metals
metals = ['Ag', 'Al', 'Au', 'Ba', 'Be', 'Bi', 'Ca', 'Cd', 'Ce', 'Co', 'Cr', 'Cs', 'Cu', 'Dy', 'Er', 'Eu', 'Fe', 'Ga', 'Gd', 'Hf', 'Hg', 'Ho', 'In', 'Ir', 'K', 'La', 'Li', 'Lu', 'Mg', 'Mn', 'Mo', 'Na', 'Nb', 'Nd', 'Ni', 'Np', 'Pb', 'Pd', 'Pr', 'Pt', 'Pu', 'Rb', 'Re', 'Rh', 'Ru', 'Sc', 'Sm', 'Sn', 'Sr', 'Tb', 'Tc', 'Th', 'Ti', 'Tl', 'Tm', 'U', 'V', 'W', 'Y', 'Yb', 'Zn', 'Zr', 'B', 'As', 'Ge', 'Sb', 'At', 'Si', 'Te']

def contains_metal(formula):
    """Check if a chemical formula contains any metal atoms"""
    for metal in metals:
        if re.search(r'(?<![a-z])' + metal + r'(?![a-z])', formula):
            return True
    return False

def generalize_formula(formula):
    """Replace specific metals with 'M' to create generalized patterns"""
    generalized = formula
    # Sort metals by length (longest first) to avoid partial replacements
    sorted_metals = sorted(metals, key=len, reverse=True)
    for metal in sorted_metals:
        # Use negative lookbehind/lookahead to ensure we don't match partial strings
        generalized = re.sub(r'(?<![a-z])' + metal + r'(?![a-z])', 'M', generalized)
    return generalized

def to_hill_subscript(prod_str):
    """Convert formula to subscript notation, handling 'M' for metals"""
    return re.sub(r"(\d+)", r"$_{\1}$", prod_str)

# Load JSON file
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)

# Build counter of metal-containing decomposition products
counter = Counter()
for entry in data.values():
    decomp = entry.get("decomposition_products", {})
    for product in decomp.keys():
        if contains_metal(product):
            generalized = generalize_formula(product)
            counter[generalized] += 1

# Get top 10 most common metal-containing products
top10 = counter.most_common(10)
products, counts = zip(*top10)

# Create DataFrame
df = pd.DataFrame({
    "product": products,
    "count": counts
})

df["pretty"] = df["product"].map(to_hill_subscript)

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
sns.barplot(
    data=df,
    x="count",
    y="product",
    palette="viridis"
)
ax.tick_params(labelsize=20)
ax.minorticks_on()
ax.tick_params(
    which='major', direction='in', length=26, width=1.8
)
ax.tick_params(
    which='minor', direction='in', length=14, width=1.6
)
for spine in ax.spines.values():
    spine.set_linewidth(2.5)
ax.tick_params(
    axis='y',
    which='both',
    length=0
)

plt.xlabel("Number of Occurrences", fontsize=22)
plt.ylabel("Decomposition Product", fontsize=22)

ax.set_yticklabels(df["pretty"], fontsize=20)

plt.tight_layout()
plt.savefig("FigureS22", dpi=1000)
plt.show()
