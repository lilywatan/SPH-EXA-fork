import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import re
import seaborn as sns
from math import pi 
from scipy.spatial import cKDTree
from matplotlib.colors import Normalize

f5 = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_1J_5000.hdf5'
# extract data from hdf5 file with optional stride # returns arrays of selected keys (based on stride), density, pressure, x, y and z coords of particles & star
def read_hdf5_data(file, stride=1):
    x_pos = []
    y_pos = []
    star_x = []
    star_y = []
    star_m = []

    with h5py.File(file, 'r') as f: 
    # extract all star mass and sound speed data
            x_pos.append(np.array(f['Step#0']['x']))
            y_pos.append(np.array(f['Step#0']['y']))
            star_x.append(np.array(f['Step#0'].attrs['star::x']))
            star_y.append(np.array(f['Step#0'].attrs['star::y']))
            star_m.append(np.array(f['Step#0'].attrs['star::m']))

    return np.array(x_pos), np.array(y_pos), star_x, star_y, star_m
# find particle with largest semi-major axis and return its index
def largest_a(x, y, star_x, star_y):
    a = np.sqrt((x - star_x)**2 + (y - star_y)**2)
    return np.argmax(a)
# calculate orbital period of particle with largest semi-major axis 
def orbital_period(x, y, star_x, star_y, star_m):
# find particle with largest semi-major axis
    i = largest_a(x, y, star_x, star_y)
    print(x.shape, y.shape)
# calculate orbital period
    a = np.sqrt((x[0][i] - star_x)**2 + (y[0][i] - star_y)**2)
    T = 2 * pi * np.sqrt(a**3 / star_m)
    print(f'Orbital period: {T}')

if __name__ == '__main__':
    x_pos, y_pos, star_x, star_y, star_m = read_hdf5_data(f5)
    orbital_period(x_pos, y_pos, star_x, star_y, star_m)