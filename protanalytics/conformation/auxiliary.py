import os
import numpy as np 
import scipy as sp 


def check_directories(directories):
  for directory in directories:
    if not os.path.exists(directory):
      os.makedirs(directory)


def angle_distance_metric(a, b):
  d = abs(a - b) % 360
  return 360 - d if d > 180 else d
  
