import sys, os, json, datetime
import numpy as np
import pandas as pd
from sklearn.neighbors import KernelDensity
from sklearn import pipeline, svm
from sklearn.cluster import DBSCAN, AffinityPropagation, KMeans
from auxiliary import *


def read_conformation_data(residue_name):
  try:
    return pd.read_csv(root_dir + filtered_dir + residue_name + ".csv")
  except: 
    return pd.read_csv(root_dir + conformations_dir + residue_name + ".csv")


def base_analytics(all_data, z_matrix, grouping_parameters):
  if simple_describe:
    all_data.describe().to_csv(logs_dir + res + "--description.csv")
  if angle_correlations:
    z_matrix.corr().to_csv(logs_dir + res + "--angle-correlations.csv")
  if parameter_sample_counts:
    samples_per_group = pd.DataFrame()
    for parameter in grouping_parameters:
      parameter = tuple(parameter)
      samples_per_group[parameter] = pd.Series({k: len(v) for k, v in all_data.groupby(parameter).groups.items()}).describe()
    samples_per_group.to_csv(logs_dir + res + "--parameter-sample-counts.csv")


def make_dirs_for_grouping(results_dir, grouping_parameters):
  for parameter in grouping_parameters:
    if not os.path.exists(results_dir + "-".join(parameter) + "/"):
      os.makedirs(results_dir + "-".join(parameter) + "/")
    

globals().update(json.load(open("settings.json")))  
log = open(logs_dir + ('{0:%Y-%m-%d %Hh %Mm}'.format(datetime.datetime.now())) + ".log", "w")

check_directories([logs_dir, results_dir, stats_dir, plots_dir])
make_dirs_for_grouping(results_dir, grouping_parameters)
overview_stats = open(stats_dir + "overview.csv", "w")


overview_stats.write("residue,kde_bandwidth,dpgmm_epsilon,dpgmm_cluster_count,affinity_cluster_count,kmeans_cluster_count\n")

for res in residues:
  print("Current working residue : " + res)
  overview_stats.write(res + ',')


  try:
    all_data = read_conformation_data(res)
  except:
    log.write(res + " has no conformation file.\n")
    continue

  with open(logs_dir + res + ".log", "w") as res_log:
    z_matrix = all_data.iloc[:,6:]

    base_analytics(all_data, z_matrix, grouping_parameters)

    res_log.write("\n\n# Clustering and Analysis \n")

    kde_bandwidth = (5 * np.mean(list(all_data.iloc[:,6:].std()))) / (len(all_data.iloc[:,6:])**(.2))
    res_log.write("KDE bandwidth : " + str(kde_bandwidth) + "\n")
    kde = KernelDensity(kernel="gaussian", bandwidth=kde_bandwidth)
    full_kde = KernelDensity(kernel="gaussian", bandwidth=kde_bandwidth)
   # dpgmm = DPGMM(n_components = 10, covariance_type = "diag")
    dbscan_epsilon = dbscan_modifier * np.sqrt(all_data.iloc[:,6:].max(axis=0).apply(lambda x : np.power(x,2)).sum(axis=0))
    res_log.write("DBSCAN epsilon : " + str(dbscan_epsilon) + "\n")
    dbscan = DBSCAN(eps=dbscan_epsilon, min_samples=5, metric='euclidean', algorithm="auto", leaf_size=30, p=None, random_state=None)
    affinity = AffinityPropagation(damping=0.8, max_iter=200, convergence_iter=15, copy=True, preference=None, affinity='euclidean', verbose=False)
    kmeans = KMeans(n_clusters=10, init='k-means++', n_init=10, max_iter=300, tol=0.0001, precompute_distances=True, verbose=0, random_state=None, copy_x=True, n_jobs=1)
    


    for parameter in grouping_parameters:
      print("Beginning parametereter: " + str(parameter))
  #    scores_file = open(results_dir + "-".join(parameter) + "/" + res + "-Scores.csv", "w")
  #    scores_file.write(",".join(list(all_data.iloc[:,:5].keys())) + ",Full KDE,Affinity Cluster,KMeans Cluster," + ",".join(list(all_data.iloc[:,5:].keys())) + "\n")

      with open(results_dir + "-".join(parameter) + "/" + res + "-parametereters.dat", "w") as parametereters_file:

        grouped = all_data.groupby(parameter, axis=0, sort=True)

        results = []
        for key, frame in grouped:
          if len(frame) < 10:
            print("Not enough samples in " + str(key))
            continue
          
          frame.index = [i for i in range(len(frame))]

          kmeans = KMeans(n_clusters=10, init='k-means++', n_init=10, max_iter=300, tol=0.0001, precompute_distances=True, verbose=0, random_state=None, copy_x=True, n_jobs=1)
          
          dbscan_epsilon = dbscan_modifier * np.sqrt(frame.iloc[:,6:].max(axis=0).apply(lambda x : np.power(x,2)).sum(axis=0))
          min_samples_for_core = len(frame)/100 if len(frame)/10 > 5 else 5
          dbscan = DBSCAN(eps=dbscan_epsilon, min_samples= min_samples_for_core, metric="euclidean", algorithm="auto", leaf_size=30, p=None, random_state=None)
          
          dbscan.fit(frame.iloc[:,6:].values)
          print("Finished " + str(key) + "  -- with n samples : " + str(len(frame)) + "\tNumber Of Clusters : " + str(len(np.unique(dbscan.labels_))) + "\tProportion Core Samples : " + str(len(dbscan.core_sample_indices_)/float(len(frame))))
    #      kmeans.fit_predict(frame.iloc[:,5:].values)
          kde.fit(frame.iloc[:,6:])
    #      dpgmm.fit(frame.iloc[:,6:])
          frame.insert(6, "Cluster KDE Score", kde.score_samples(frame.iloc[:,6:].values))
    #      frame.insert(6, "KMean Cluster", kmeans.labels_)
          frame.insert(6, "DBSCAN Cluster", dbscan.labels_)

          # frame.to_csv(results_dir + "-".join(parameter) + "/" + res + "-scores--" + "-".join(map(str, key)) + ".csv", index = False)
          results.append(frame)
    #      parametereters_file.write(str(key) + ",")
    #      parametereters_file.write(str(kde.metric) + ",")
    #      parametereters_file.write(str(dpgmm.weights_) + "\n")

        full_result = pd.concat(results, axis = 0, join = "outer", ignore_index = True)
        full_result.to_csv(results_dir + "-".join(parameter) + "/" + res + "--" + "-".join(parameter) + "-scores--full.csv", index = False)

