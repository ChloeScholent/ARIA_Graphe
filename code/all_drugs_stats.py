import pandas as pd
import networkx as nx
from itertools import combinations

# ============================================================
# --- Load CSV Files ---
# ============================================================
drug_protein_df = pd.read_csv("csv/bio-decagon-targets-all.csv")  # columns: STITCH, Gene
drug_drug_df = pd.read_csv("csv/bio-decagon-combo.csv")          # columns: STITCH 1, STITCH 2, Polypharmacy Side Effect
ppi_df = pd.read_csv("csv/bio-decagon-ppi.csv")                  # columns: Gene 1, Gene 2, interaction_score

# ============================================================
# --- Extract Unique Drugs and Proteins ---
# ============================================================
drugs = drug_protein_df['STITCH'].unique()
proteins = drug_protein_df['Gene'].unique()
drug_set = set(drugs)

# ============================================================
# --- Filter Relevant Drug–Drug Combos ---
# ============================================================
filtered_dd = drug_drug_df[
    drug_drug_df['STITCH 1'].isin(drug_set) &
    drug_drug_df['STITCH 2'].isin(drug_set)
]

# ============================================================
# --- Map Each Drug → Its Protein Targets ---
# ============================================================
drug_to_proteins = {
    drug: set(drug_protein_df.loc[drug_protein_df['STITCH'] == drug, 'Gene'])
    for drug in drugs
}

# ============================================================
# --- Prepare All Drug Pairs ---
# ============================================================
all_pairs = set(tuple(sorted(pair)) for pair in combinations(drugs, 2))
side_effect_pairs = set(tuple(sorted((row['STITCH 1'], row['STITCH 2']))) for _, row in filtered_dd.iterrows())
no_side_effect_pairs = all_pairs - side_effect_pairs

print(f"Total drugs: {len(drugs)}")
print(f"Total drug pairs: {len(all_pairs)}")
print(f"Side-effect pairs: {len(side_effect_pairs)}")
print(f"No-side-effect pairs: {len(no_side_effect_pairs)}")

# ============================================================
# --- Protein–Protein Interactions Set for Quick Lookup ---
# ============================================================
ppi_edges = {(min(u, v), max(u, v)) for u, v in zip(ppi_df['Gene 1'], ppi_df['Gene 2'])}

# ============================================================
# --- Compute Metrics ---
# ============================================================
def compute_metrics(pairs):
    shared_list = []
    ppi_list = []
    for d1, d2 in pairs:
        proteins_u = drug_to_proteins.get(d1, set())
        proteins_v = drug_to_proteins.get(d2, set())
        shared = proteins_u & proteins_v
        shared_list.append(len(shared))
        ppi_count = sum((min(p1, p2), max(p1, p2)) in ppi_edges for p1 in proteins_u for p2 in proteins_v)
        ppi_list.append(ppi_count)
    return shared_list, ppi_list

# Side-effect combos
shared_side, ppi_side = compute_metrics(side_effect_pairs)
# No-side-effect combos
shared_none, ppi_none = compute_metrics(no_side_effect_pairs)

# ============================================================
# --- Compute Averages ---
# ============================================================
avg_shared_side = sum(shared_side)/len(shared_side) if shared_side else 0
avg_ppi_side = sum(ppi_side)/len(ppi_side) if ppi_side else 0

avg_shared_none = sum(shared_none)/len(shared_none) if shared_none else 0
avg_ppi_none = sum(ppi_none)/len(ppi_none) if ppi_none else 0

print("\n=== Protein Interaction & Shared Target Stats (All Drugs) ===")
print(f"Average shared proteins (side-effect combos): {avg_shared_side:.2f}")
print(f"Average PPI interactions (side-effect combos): {avg_ppi_side:.2f}")
print(f"Average shared proteins (no-side-effect combos): {avg_shared_none:.2f}")
print(f"Average PPI interactions (no-side-effect combos): {avg_ppi_none:.2f}")
