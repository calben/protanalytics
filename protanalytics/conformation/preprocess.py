import sys, os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from auxiliary import *

globals().update(json.load(open("settings.json")))

def read_conformation_data(residue_name):
  return pd.read_csv(root_dir + conformations_dir + residue_name + ".csv")


def filter_conformations():
  for res in residues:
    print("Current working residue : " + res)

    try:
      all_data = read_conformation_data(res)
    except:
      print(res + " has no conformation file.\n")
      continue

    metadata = all_data.ix[:,:5]
    z_matrix = all_data.iloc[:,5:]
    z_matrix = z_matrix.applymap(np.radians)
    z_sin_cos_matrix = radian_matrix_to_sin_cos_matrix(z_matrix)
    filtered_matrix = pd.concat([metadata, z_sin_cos_matrix], axis=1, join='outer', join_axes=None, ignore_index=False,
         keys=None, levels=None, names=None, verify_integrity=True)

    filtered_matrix.to_csv(filtered_dir + res + ".csv", index = True, index_label = "Label")

filter_conformations()
