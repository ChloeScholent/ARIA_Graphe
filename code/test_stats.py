import csv
import numpy as np
from pprint import pprint as print
import pandas as pd

file_path = "csv/bio-decagon-combo.csv"
file_path_2 = "csv/bio-decagon-targets-all.csv"

drug_list = []
drug = []

with open(file_path, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:

                drug_list.append(row["STITCH 1"])
                drug_list.append(row["STITCH 2"])



drug_dict = list(dict.fromkeys(drug_list))
print(len(drug_dict))

print(len(list(set(drug_list))))


with open(file_path_2, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:

                drug.append(row["STITCH"])




drug_dict_2 = list(dict.fromkeys(drug))
print(len(drug_dict_2))

print(len(list(set(drug))))