import pandas as pd
import networkx as nx
import os
from itertools import combinations

# ============================================================
# --- File Paths ---
# ============================================================
drug_dir = "csv/random_drug_sets_80"     # Folder with your 15 random_80_drugs_set_X.csv
drug_drug_file = "csv/bio-decagon-combo.csv"
ppi_file = "csv/bio-decagon-ppi.csv"
output_summary = "csv/drug_set_stats_summary.csv"

# ============================================================
# --- Load Shared Data ---
# ============================================================
print("Loading large CSVs...")

drug_drug_df = pd.read_csv(drug_drug_file)  # STITCH 1, STITCH 2, Polypharmacy Side Effect
ppi_df = pd.read_csv(ppi_file)              # Gene 1, Gene 2, interaction_score

# Convert PPI into a quick lookup set
ppi_edges = {(min(u, v), max(u, v)) for u, v in zip(ppi_df['Gene 1'], ppi_df['Gene 2'])}

# ============================================================
# --- Prepare Output Storage ---
# ============================================================
results_summary = []

# ============================================================
# --- Iterate Over Each Drug Subset File ---
# ============================================================
for file in sorted(os.listdir(drug_dir)):
    if not file.endswith(".csv"):
        continue

    path = os.path.join(drug_dir, file)
    print(f"\nProcessing {file}...")

    # Load subset
    drug_protein_df = pd.read_csv(path)

    # Build sets
    drugs = drug_protein_df['STITCH'].unique()
    proteins = drug_protein_df['Gene'].unique()
    drug_set = set(drugs)

    # Filter drug–drug combos relevant to this set
    filtered_dd = drug_drug_df[
        drug_drug_df['STITCH 1'].isin(drug_set) &
        drug_drug_df['STITCH 2'].isin(drug_set)
    ]

    # Map each drug → its protein targets
    drug_to_proteins = {
        drug: set(drug_protein_df.loc[drug_protein_df['STITCH'] == drug, 'Gene'])
        for drug in drugs
    }

    # Compute all possible pairs
    all_pairs = set(tuple(sorted(pair)) for pair in combinations(drugs, 2))

    # Known side-effect pairs
    side_effect_pairs = set(tuple(sorted((row['STITCH 1'], row['STITCH 2'])))
                            for _, row in filtered_dd.iterrows())

    # No-side-effect pairs
    no_side_effect_pairs = all_pairs - side_effect_pairs

    # Stats containers
    side_effect_stats = []
    no_side_effect_stats = []

    # --- Compute for side-effect pairs ---
    for d1, d2 in side_effect_pairs:
        proteins_u = drug_to_proteins.get(d1, set())
        proteins_v = drug_to_proteins.get(d2, set())
        shared_proteins = proteins_u & proteins_v
        num_shared = len(shared_proteins)
        # PPI count between targets
        num_ppi_links = sum(
            (min(p1, p2), max(p1, p2)) in ppi_edges
            for p1 in proteins_u for p2 in proteins_v
        )
        side_effect_stats.append((num_shared, num_ppi_links))

    # --- Compute for no-side-effect pairs ---
    for d1, d2 in no_side_effect_pairs:
        proteins_u = drug_to_proteins.get(d1, set())
        proteins_v = drug_to_proteins.get(d2, set())
        shared_proteins = proteins_u & proteins_v
        num_shared = len(shared_proteins)
        num_ppi_links = sum(
            (min(p1, p2), max(p1, p2)) in ppi_edges
            for p1 in proteins_u for p2 in proteins_v
        )
        no_side_effect_stats.append((num_shared, num_ppi_links))

    # Compute averages
    avg_shared_side = (sum(s for s, _ in side_effect_stats) / len(side_effect_stats)) if side_effect_stats else 0
    avg_ppi_side = (sum(p for _, p in side_effect_stats) / len(side_effect_stats)) if side_effect_stats else 0

    avg_shared_none = (sum(s for s, _ in no_side_effect_stats) / len(no_side_effect_stats)) if no_side_effect_stats else 0
    avg_ppi_none = (sum(p for _, p in no_side_effect_stats) / len(no_side_effect_stats)) if no_side_effect_stats else 0

    # Store results
    results_summary.append({
        "file": file,
        "num_drugs": len(drugs),
        "num_side_effect_pairs": len(side_effect_pairs),
        "num_no_side_effect_pairs": len(no_side_effect_pairs),
        "avg_shared_proteins_side_effect": avg_shared_side,
        "avg_ppi_side_effect": avg_ppi_side,
        "avg_shared_proteins_no_side_effect": avg_shared_none,
        "avg_ppi_no_side_effect": avg_ppi_none
    })

# ============================================================
# --- Save Summary ---
# ============================================================
summary_df = pd.DataFrame(results_summary)
summary_df.to_csv(output_summary, index=False)

print(f"\n✅ Summary saved to: {output_summary}")
print(summary_df)
