import pandas as pd

# Load the CSV
df = pd.read_csv("csv/bio-decagon-combo.csv")

# Step 1: Normalize side effect names (lowercase, strip spaces)
df["Side Effect Name"] = df["Side Effect Name"].str.lower().str.strip()

# Step 2: Find all pairs with "hypertension"
pairs_with_hypertension = df[df["Side Effect Name"] == "pulmonary hypertension"]

# Count unique drug pairs with hypertension
num_hypertension = pairs_with_hypertension[["STITCH 1", "STITCH 2"]].drop_duplicates().shape[0]
print(f"Number of pairs with hypertension: {num_hypertension}")

# Step 3: Among those pairs, find how many also have "fever"
# Create sets of pairs for each effect
hypertension_pairs = set(zip(pairs_with_hypertension["STITCH 1"], pairs_with_hypertension["STITCH 2"]))

fever_rows = df[df["Side Effect Name"] == "hay fever"]
fever_pairs = set(zip(fever_rows["STITCH 1"], fever_rows["STITCH 2"]))

# Intersection: pairs that have both hypertension and fever
both_effects_pairs = hypertension_pairs.intersection(fever_pairs)

print(f"Number of pairs with both hypertension and fever: {len(both_effects_pairs)}")

print(f'{len(both_effects_pairs)*100/num_hypertension:.2f}')