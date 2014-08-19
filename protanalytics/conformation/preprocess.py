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
    columns_to_keep = map(lambda x : True if x > angle_minimum_std else False, z_matrix.apply(np.std).values)

    z_sin_matrix = z_matrix.apply(np.sin)
    z_sin_matrix.columns = map(lambda x : x + "--SIN", z_matrix.columns)
    z_cos_matrix = z_matrix.apply(np.cos)
    z_cos_matrix.columns = map(lambda x : x + "--COS", z_matrix.columns)

    filtered_matrix = pd.concat([metadata, z_sin_matrix, z_cos_matrix], axis=1, join='outer', join_axes=None, ignore_index=False,
         keys=None, levels=None, names=None, verify_integrity=True)

    filtered_matrix.to_csv(filtered_dir + res + ".csv", index = True, index_label = "Label")

filter_conformations()
