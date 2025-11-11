import csv
import numpy as np

# === File paths ===
file_path = "csv/bio-decagon-combo.csv"
file_path_protein = "csv/bio-decagon-targets-all.csv"

# === UMLS codes for target side effects ===
FEVER_CODE = "C0018621"
HYPERTENSION_CODE = "C0020542"

# === Data structures ===
pair_to_side_effects = {}
drug_pair = []
test = []

# === STEP 1: Read all drug-side effect combinations ===
with open(file_path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        pair = tuple(sorted([row["STITCH 1"], row["STITCH 2"]]))  # order-independent pair
        side_effect = row["Polypharmacy Side Effect"]

        if pair not in pair_to_side_effects:
            pair_to_side_effects[pair] = set()
        pair_to_side_effects[pair].add(side_effect)

# === STEP 2: Select pairs that have BOTH fever & hypertension ===
selected_pairs = [
    pair for pair, effects in pair_to_side_effects.items()
    if FEVER_CODE in effects and HYPERTENSION_CODE in effects
]

print(f"✅ Number of drug pairs with both fever and hypertension: {len(selected_pairs)}")

# === STEP 3: Extract unique drugs involved in those pairs ===
drug_list = sorted(set([drug for pair in selected_pairs for drug in pair]))

print(f"✅ Number of unique drugs involved: {len(drug_list)}")

# === STEP 4: Load protein associations ===
drug_combo = np.loadtxt(
    file_path_protein,
    delimiter=",",
    dtype=[('col1', '<U26'), ('col2', 'i4')],
    skiprows=1
)

# === STEP 5: Write matching drugs and their proteins ===
output_path = "csv/fever_hyperT_drug_protein.csv"

with open(output_path, "w", newline='') as output_file:
    writer = csv.writer(output_file)
    writer.writerow(["Drug", "Protein"])  # optional header

    for row in drug_combo:
        if row[0] in drug_list:
            writer.writerow([row[0], row[1]])
            test.append(row[0])

# === STEP 6: Deduplicate and report ===
test = list(dict.fromkeys(test))
print(f"✅ Unique drugs written to {output_path}: {len(test)}")
