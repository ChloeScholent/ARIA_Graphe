import csv
import numpy as np
from pprint import pprint as print
import pandas as pd

file_path = "csv/bio-decagon-combo.csv"
file_path_protein = "csv/bio-decagon-targets-all.csv"
diff_drug_path = "csv/bio-decagon-combo.csv"

drug_list = []
drug_pair = []
diff_drug_pair = []
diff_drug_list = []
test = []
test_2 = []


####################################
#####SIMILAR SIDE EFFECT COMBO#####
##################################

with open(file_path, newline='') as f:
    reader = csv.DictReader(f)
    x = 0
    for row in reader:
            if row["Polypharmacy Side Effect"] == "C0006826" and x<30:
                x+=1
                drug_list.append(row["STITCH 1"])
                drug_list.append(row["STITCH 2"])
                pair = {row["STITCH 1"], row["STITCH 2"]}
                drug_pair.append(pair)

drug_list = list(dict.fromkeys(drug_list))

# print(len(drug_list))
# print(drug_pair)
# print(len(drug_pair))

drug_combo = np.loadtxt(file_path_protein, delimiter=",", dtype=[('col1', '<U26'), ('col2', 'i4')], skiprows= 1)

with open("csv/similar_drug_protein.csv", "w", newline='') as output_file:
    writer = csv.writer(output_file)

    for row in drug_combo:
        for drug in drug_list:
            if row[0] == drug:
                writer.writerow([row[0], row[1]])
                test.append(row[0])
        
test = list(dict.fromkeys(test))
# print(len(test))



######################################
#####DIFFERENT SIDE EFFECT COMBO#####
####################################

diff_drug_combo = np.loadtxt(diff_drug_path, delimiter=",", dtype=[('col1', '<U26'), ('col2', '<U26'),('col3', '<U26'),('col4', '<U26'),], skiprows= 1)

for row in diff_drug_combo:
    # Convert to Python strings using str()
    drug1 = str(row[0])
    drug2 = str(row[1])
    diff_drug_list.append(drug1)
    diff_drug_list.append(drug2)
    pair = {drug1, drug2}
    diff_drug_pair.append(pair)

diff_drug_list = list(dict.fromkeys(diff_drug_list))

print(diff_drug_pair)
# print(len(diff_drug_list))
# print(diff_drug_list)

with open("csv/different_drug_protein.csv", "w", newline='') as output_file:
    writer = csv.writer(output_file)

    for row in drug_combo:
        for drug in diff_drug_list:
            if row[0] == drug:
                writer.writerow([row[0], row[1]])
                test_2.append(row[0])

test_2 = list(dict.fromkeys(test_2))
# print(len(test_2))






