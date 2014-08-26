import os
import numpy as np 
import scipy as sp 
import pandas as pd


def check_directories(directories):
  for directory in directories:
    if not os.path.exists(directory):
      os.makedirs(directory)


def angle_distance_metric(a, b):
  d = abs(a - b) % 360
  return 360 - d if d > 180 else d
  

def z_angle_matrix_to_z_sin_cos_matrix(z_angle_matrix):
  z_sin_matrix = z_angle_matrix.apply(np.sin)
  z_sin_matrix.columns = map(lambda x : x + "--SIN", z_angle_matrix.columns)
  z_cos_matrix = z_angle_matrix.apply(np.cos)
  z_cos_matrix.columns = map(lambda x : x + "--COS", z_angle_matrix.columns)
  return pd.concat([z_sin_matrix, z_cos_matrix], axis=1, join='outer', join_axes=None, ignore_index=False, keys=None, levels=None, names=None, verify_integrity=True)

def z_sin_cos_matrix_to_z_angle_matrix(z_sin_cos_matrix):
  z_sin_matrix = z_sin_cos_matrix[0:(len(z_sin_cos_matrix)/2)]
  z_cos_matrix = z_sin_cos_matrix[(len(z_sin_cos_matrix)/2):0]
