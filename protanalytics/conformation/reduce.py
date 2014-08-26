import pandas as pd 
import numpy as np 
import json, sys
from auxiliary import *


globals().update(json.load(open("settings.json")))


def read_results_data(residue_name, parameter):
  print(results_dir + "-".join(parameter) + "/" + res + "-scores--full.csv")
  return pd.read_csv(results_dir + "-".join(parameter) + "/" + res + "-scores--full.csv")


def analyse_group(df):
  z_matrix = df.iloc[:,8:]
  means = z_matrix.mean(axis=0, skipna=False, level=None, numeric_only=None)
  medians = z_matrix.median(axis=0, skipna=False, level=None, numeric_only=None)
  kurtoses = z_matrix.kurt(axis=0, numeric_only=None)
  stds = z_matrix.std(axis=0, skipna=False, level=None, numeric_only=None)
  df = df.sort("Cluster KDE Score", ascending=False)
  best = df.iloc[0]
  return np.concatenate([best.values, means.values, medians.values, kurtoses.values, stds.values, [len(df)]])
  

for res in residues:
  for param in grouping_parameters:
    print("Current working residue : " + res + " " + str(param))
    try:
      all_data = read_results_data(res, param)
    except:
      print(res + " has no results file.\n")
      continue

    metadata = all_data.iloc[:,:8]
    sin_cos_matrix = all_data.iloc[:,8:]
    radian_matrix = sin_cos_matrix_to_radian_matrix(sin_cos_matrix)
    angle_matrix = radian_matrix.applymap(np.degrees)
    all_data = pd.concat([metadata, angle_matrix], axis=1, join='outer', join_axes=None, ignore_index=False, keys=None, levels=None, names=None, verify_integrity=True)

    output = open(results_dir + "-".join(param) + "/" + res + "-reduced.csv", "w")
    output.write(",".join(all_data.columns) + ",")
    output.write(",".join(map(lambda x : x + "-mean", all_data.columns[8:])) + ",")
    output.write(",".join(map(lambda x : x + "-median", all_data.columns[8:])) + ",")
    output.write(",".join(map(lambda x : x + "-kurt", all_data.columns[8:])) + ",")
    output.write(",".join(map(lambda x : x + "-std", all_data.columns[8:])) + ",Samples In Cluster\n")

#    to get the necessary stats    
#    grouped.get_group((0,0,0,0)).describe().iloc[1:4,8:]

    grouping_param_with_cluster = param[:]
    grouping_param_with_cluster.append("DBSCAN Cluster")
    grouped = all_data.groupby(grouping_param_with_cluster)

    for k, v in grouped:
      output.write(",".join(map(str, analyse_group(v))) + "\n")
