import pandas as pd 
import numpy as np 
import json, sys

globals().update(json.load(open("settings.json")))


def read_results_data(residue_name, parameter):
  print(results_dir + "-".join(parameter) + "/" + res + "-scores--full.csv")
  return pd.read_csv(results_dir + "-".join(parameter) + "/" + res + "-scores--full.csv")


def analyse_group(df):
  z_matrix = df.iloc[:,8:]
  means = z_matrix.mean(axis=0)
  stds = z_matrix.std(axis=0)
  df = df.sort("Cluster KDE Score", ascending=False)
  best = df.iloc[0]
  return np.concatenate([best.values, means.values, stds.values, [len(df)]])
  

for res in residues:
  for param in grouping_parameters:
    print("Current working residue : " + res + " " + str(param))
    try:
      all_data = read_conformation_data(res, param)
    except:
      print(res + " has no results file.\n")
      continue

    output = open(results_dir + "-".join(param) + "/" + res + "-reduced.csv", "w")
    output.write(",".join(all_data.columns) + ",")
    output.write(",".join(map(lambda x : x + "-mean", all_data.columns[8:])) + ",")
    output.write(",".join(map(lambda x : x + "-std", all_data.columns[8:])) + ",Samples In Cluster\n")

#    to get the necessary stats    
#    grouped.get_group((0,0,0,0)).describe().iloc[1:4,8:]

    grouping_param_with_cluster = param[:]
    grouping_param_with_cluster.append("DBSCAN Cluster")
    grouped = all_data.groupby(grouping_param_with_cluster)

    for k, v in grouped:
      output.write(",".join(map(str, analyse_group(v))) + "\n")
