from __future__ import annotations

import gzip
import json
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

# Load data
with gzip.open("All_qmof_results.json.gz", "rt") as f:
    data = json.load(f)


# Function to extract pattern from mofid
def extract_pattern(mofid):
    if mofid is None or mofid == "null" or not isinstance(mofid, str):
        return None

    if " MOFid" not in mofid:
        return None

    before_mofid = mofid.split(" MOFid")[0]

    if "." not in before_mofid:
        return None

    parts = before_mofid.split(".", 1)
    if len(parts) != 2:
        return None

    node = parts[0]
    linker = parts[1]

    return (node, linker)


# Function to convert frac_composition to comparable format
def frac_comp_to_tuple(frac_comp):
    if frac_comp is None:
        return None
    return tuple(sorted(frac_comp.items()))


# Group entries by pattern
pattern_to_entries = defaultdict(list)

for qmof_id, entry in data.items():
    mofid = entry.get("mofid")
    pattern = extract_pattern(mofid)

    if pattern is not None:
        pattern_to_entries[pattern].append(
            {
                "qmof_id": qmof_id,
                "mofid": mofid,
                "ehull": entry.get("ehull"),
                "synthesizable": entry.get("synthesizable"),
                "frac_composition": entry.get("frac_composition"),
            }
        )

# Find and print matches (patterns with more than 1 entry AND same frac_composition)
print("MOFid Polymorph Matches (same node/linker AND same composition):")
print("=" * 100)

polymorph_count = 0
polymorph_groups = 0


for (node, linker), entries in sorted(pattern_to_entries.items()):
    if len(entries) > 1:  # Multiple entries with same pattern
        # Check if all have the same frac_composition
        frac_comps = [frac_comp_to_tuple(e["frac_composition"]) for e in entries]
        # Skip if any frac_composition is None
        if None in frac_comps:
            continue
        # Check if all frac_compositions are identical
        if len(set(frac_comps)) == 1:  # All the same
            polymorph_groups += 1
            polymorph_count += len(entries)
            print(f"\nPattern: Node=[{node}], Linker=[{linker}]")
            print(f"Number of polymorphs: {len(entries)}")
            # Print composition
            frac_comp = entries[0]["frac_composition"]
            comp_str = ", ".join(
                [f"{k}: {v:.3f}" for k, v in sorted(frac_comp.items())]
            )
            print(f"Composition: {comp_str}")
            print("-" * 100)
            for entry in entries:
                synth_str = (
                    "Synthesized" if entry["synthesizable"] else "Not Synthesized"
                )
                print(
                    f"  {entry['qmof_id']:<20} | {synth_str:<20} | E_hull: {entry['ehull']:<10.4f} | {entry['mofid']}"
                )

# Summary statistics
total_patterns = len(pattern_to_entries)
patterns_with_matches = sum(
    1 for entries in pattern_to_entries.values() if len(entries) > 1
)
total_matched_entries = sum(
    len(entries) for entries in pattern_to_entries.values() if len(entries) > 1
)

print("\n" + "=" * 100)
print(f"Total unique patterns: {total_patterns}")
print(f"Patterns with multiple entries: {patterns_with_matches}")
print(f"True polymorph groups (same composition): {polymorph_groups}")
print(f"Total polymorphs: {polymorph_count}")
print(
    f"Non-polymorphs (same MOFid pattern, different composition): {total_matched_entries - polymorph_count}"
)


# Collect ehull ranges for each polymorph group
ehull_ranges = []

for entries in pattern_to_entries.values():
    if len(entries) > 1:
        frac_comps = [frac_comp_to_tuple(e["frac_composition"]) for e in entries]

        if None in frac_comps:
            continue

        if len(set(frac_comps)) == 1:  # True polymorphs
            # Get all ehull values for this group
            ehulls = [e["ehull"] for e in entries if e["ehull"] is not None]

            if len(ehulls) > 1:
                ehull_range = max(ehulls) - min(ehulls)
                ehull_ranges.append(ehull_range)

# Plot histogram
fig, ax = plt.subplots(figsize=(10, 6))

ax.hist(ehull_ranges, bins=50, color="#6495ED", alpha=0.7, edgecolor="black")

ax.set_xlabel(r"Δ$E_{\mathrm{hull}}$ Range (eV/atom)", fontsize=16)
ax.set_ylabel("Number of Polymorph Groups", fontsize=16)
ax.set_title(
    r"Distribution of Δ$E_{\mathrm{hull}}$ Ranges for Polymorph Groups", fontsize=16
)

ax.tick_params(which="major", direction="in", length=10, width=1.5)
ax.tick_params(which="minor", direction="in", length=5, width=1.5)
ax.tick_params(labelsize=12)
ax.minorticks_on()

for spine in ax.spines.values():
    spine.set_linewidth(1.5)

# Add statistics on plot
mean_range = np.mean(ehull_ranges)
median_range = np.median(ehull_ranges)
# ax.axvline(mean_range, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_range:.4f}')
# ax.axvline(median_range, color='green', linestyle='--', linewidth=2, label=f'Median: {median_range:.4f}')
# ax.legend(fontsize=12)

plt.tight_layout()
plt.savefig("FigureS8.png", dpi=300, bbox_inches="tight")
plt.show()

print("\nE_hull Range Statistics:")
print(f"Mean range: {mean_range:.4f} eV/atom")
print(f"Median range: {median_range:.4f} eV/atom")
print(f"Min range: {min(ehull_ranges):.4f} eV/atom")
print(f"Max range: {max(ehull_ranges):.4f} eV/atom")
