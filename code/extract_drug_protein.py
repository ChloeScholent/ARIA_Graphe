import csv
import numpy as np

csv_file_1 = "csv/bio-decagon-targets-all.csv"
csv_file_2 = "csv/bio-decagon-ppi.csv"

protein_list = []
drug_list = []

# with open(csv_file_1, newline='') as f:
#     reader = csv.DictReader(f)
#     for row in reader:
#         drug_list.append(row['STITCH'])
#         protein_list.append(row['Gene'])

# short_drug_list = list(dict.fromkeys(drug_list))
# print(len(short_drug_list))

# short_protein_list = list(dict.fromkeys(protein_list))
# print(len(short_protein_list))

# with open("short_protein.csv", "a", newline='') as file:
#     writer = csv.writer(file)
#     with open(csv_file_2, newline='') as f:
#         reader = csv.DictReader(f)
#         for row in reader:
#             for protein in short_protein_list:
#                 if protein == row['Gene 1']:
#                     writer.writerow([row['Gene 1'], row['Gene 2']])  # appends a single new row


my_data= np.genfromtxt(csv_file_2, delimiter=",")


