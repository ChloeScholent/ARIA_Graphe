import pandas as pd


df = pd.read_csv("csv/bio-decagon-combo.csv")

df["Side Effect Name"] = df["Side Effect Name"].str.lower().str.strip()

pairs_with_hypertension = df[df["Side Effect Name"] == "pulmonary hypertension"]

num_hypertension = pairs_with_hypertension[["STITCH 1", "STITCH 2"]].drop_duplicates().shape[0]
print(f"Number of pairs with hypertension: {num_hypertension}")

hypertension_pairs = set(zip(pairs_with_hypertension["STITCH 1"], pairs_with_hypertension["STITCH 2"]))

fever_rows = df[df["Side Effect Name"] == "hay fever"]
fever_pairs = set(zip(fever_rows["STITCH 1"], fever_rows["STITCH 2"]))

both_effects_pairs = hypertension_pairs.intersection(fever_pairs)

print(f"Number of pairs with both hypertension and fever: {len(both_effects_pairs)}")

print(f'{len(both_effects_pairs)*100/num_hypertension:.2f}')