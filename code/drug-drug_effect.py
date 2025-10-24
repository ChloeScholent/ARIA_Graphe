import csv
import numpy as np


file_path = "csv/bio-decagon-combo.csv"
file_path_protein = "csv/bio-decagon-targets.csv"
x = 0
drug_list = []
test = []

with open(file_path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
            if row["Side Effect Name"] == "eruption" and x<30:
                x+=1
                drug_list.append(row["STITCH 1"])
                drug_list.append(row["STITCH 2"])

drug_list = list(dict.fromkeys(drug_list))

drug_combo = np.loadtxt(file_path_protein, delimiter=",", dtype=[('col1', '<U26'), ('col2', 'i4')], skiprows= 1)

with open("csv/similar_drug_protein.csv", "w", newline='') as output_file:
    writer = csv.writer(output_file)

    for row in drug_combo:
        if row[0] in drug_list:
            writer.writerow([row[0], row[1]])
            test.append(row[0])
    
test = list(dict.fromkeys(test))
print(len(test))
print(len(drug_list))
4543