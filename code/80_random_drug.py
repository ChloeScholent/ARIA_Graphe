import pandas as pd
import random
import os

# --- Load your master CSV ---
df = pd.read_csv("csv/bio-decagon-targets-all.csv")  # must contain 'STITCH' column

# --- Get all unique drugs ---
unique_drugs = list(df['STITCH'].unique())

# --- Shuffle for randomness ---
random.shuffle(unique_drugs)

# --- Define how many drugs per subset and total subsets ---
drugs_per_file = 80
num_files = 15

# --- Safety check ---
total_needed = drugs_per_file * num_files
if total_needed > len(unique_drugs):
    raise ValueError(f"Not enough unique drugs ({len(unique_drugs)}) for {num_files} files of {drugs_per_file} each.")

# --- Create output directory ---
os.makedirs("csv/random_drug_sets_80", exist_ok=True)

# --- Split into disjoint sets ---
for i in range(num_files):
    start_idx = i * drugs_per_file
    end_idx = start_idx + drugs_per_file
    drug_subset = unique_drugs[start_idx:end_idx]
    
    # Filter dataframe for these drugs
    subset_df = df[df['STITCH'].isin(drug_subset)]
    
    # Save to new CSV
    output_path = f"csv/random_drug_sets_80/random_80_drugs_set_{i+1}.csv"
    subset_df.to_csv(output_path, index=False)
    
    print(f"✅ Saved {len(drug_subset)} unique drugs with {len(subset_df)} interactions → {output_path}")

print("\nAll 15 CSV files created successfully. Each file has 100 unique, non-overlapping drugs.")
