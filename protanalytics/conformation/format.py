import sys, os, json
from auxiliary import *


root_dir = ""

globals().update(json.load(open("settings.json")))  


labelled = dict(zip(residues, [False] * 20))

residues = dict(zip(residues, range(len(residues))))
secondary_structures = dict(zip(secondary_structures, range(len(secondary_structures))))

recognised_symbols = dict(list(residues.items()) + list(secondary_structures.items()))            
check_directories([logs_dir, results_dir, stats_dir, plots_dir, filtered_dir, conformations_dir])

print(recognised_symbols)

def convert_recognised_symbol_to_num(line, residues, secondary_structures):
  for k, v in recognised_symbols.items():
    line = line.replace(str(k), str(v))
  return line


def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False


def initialise_conformation_files(root_dir):
  conformation_files = dict()
  for k in residues:
    conformation_files[k] = (open(conformations_dir + k + ".csv", "w"))
  return conformation_files


conformation_files = initialise_conformation_files(root_dir)

with open(root_dir + "data/conformations.txt") as f:
  for line in f:
    residue_type = line[0:3]
    try: residues[residue_type] # checks if first three chars represent a residue
    except(KeyError): continue # skips comment lines and other non residue lines
    if (residue_type == "GLY" or residue_type == "PRO" or residue_type == "ALA"): continue
    
    line = convert_recognised_symbol_to_num(line, residues, secondary_structures)
    try:
        (type_to_secstruc, after_ca) = line.split("CA: ", 1)
        (z_matrix, front_to_interacting) = after_ca.split(" adj.: ", 1)
    except ValueError: continue

    z_matrix = " ".join(z_matrix.split()[7:]) # kills the leading CA and CA-CB    
    
    if not labelled[residue_type]:
        conformation_files[residue_type].write("Phi,Psi,SecStruc,Front,Back," + 
                ",".join([s[:-1] for s in z_matrix.split() if not is_number(s)]) + "\n")
        labelled[residue_type] = True
    
    z_matrix = [s for s in z_matrix.split() if is_number(s)] # z matrix only numbers from text
    z_matrix = z_matrix[::-3][::-1] # z matrix only angles

    front_to_back = front_to_interacting[0: front_to_interacting.index(" interacts with:")].strip().replace(" ", ",")
    final_line = ",".join(map(str, type_to_secstruc.split()[2:])) + "," + front_to_back + \
            "," + ",".join(map(str, z_matrix))
    if not final_line.replace(",", "").replace(" ", "").replace(".", "").replace("-", "").isdigit():
        continue
    conformation_files[residue_type].write(final_line + "\n")
