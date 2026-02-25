import json
import pandas as pd
import matplotlib.pyplot as plt
import gzip

# Load r2scan results (25 entries)
with gzip.open("r2scan_qmof_results.json.gz", "rt") as f:
    r2scan_data = json.load(f)

# Load PBE results (full database)
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    pbe_data = json.load(f)

# Build matching dataset
matched_data = []
for qmof_id, r2scan_entry in r2scan_data.items():
    if qmof_id in pbe_data:
        matched_data.append({
            'qmof_id': qmof_id,
            'pbe_eform': pbe_data[qmof_id]['formation_energy'],
            'r2scan_eform': r2scan_entry['formation_energy']
        })
    else:
        print(f"Warning: {qmof_id} not found in All_qmof_results.json")

# Convert to DataFrame
df = pd.DataFrame(matched_data)

print(f"Matched {len(df)} entries out of {len(r2scan_data)} r2scan entries")

# Create scatter plot
fig, ax = plt.subplots(figsize=(8, 7))
plt.scatter(df['pbe_eform'], df['r2scan_eform'], s=50, c='k', alpha=0.7)

# Plot y = x line
plt.plot([-1, 1], [-1, 1], linestyle='--', color='black', linewidth=2)

# Labels and styling
plt.xlabel(r'PBE-D3(BJ) Δ$E_{\mathrm{form}}$ (eV/atom)', fontsize=22)
plt.ylabel(r'r$^{\mathrm{2}}$SCAN-D4 Δ$E_{\mathrm{form}}$ (eV/atom)', fontsize=22)

ax.tick_params(which='major', direction='in', length=26, width=1.8)
ax.tick_params(which='minor', direction='in', length=14, width=1.6)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)

ax.tick_params(labelsize=20)
ax.minorticks_on()
ax.tick_params(axis='x', which='major', pad=12)
ax.tick_params(axis='y', which='major', pad=12)

plt.xlim([-0.7, 0.2])
plt.ylim([-0.7, 0.2])

ax.text(-0.2, 1.05, "A", transform=ax.transAxes, fontsize=24,
        fontweight="bold", va="top", ha="left")

plt.tight_layout()
plt.savefig("FigureS5A", dpi=300)
plt.show()
