import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import re

f5 = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_1J_5000.hdf5'
f5_double = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_2J_5000.hdf5'
f5_half = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_05J_5000.hdf5'
output = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/5000.txt'
plots = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/plots/'

# constants
G = 1 
G2 = G * G

# extract data from hdf5 file with optional stride 
def read_hdf5_data(file, stride=1):
    densities = []
    pressures = []
    x_pos = []
    y_pos = []
    z_pos = []
    times = []
    star_x = []
    star_y = []
    star_z = []

    with h5py.File(file, 'r') as f: 
        step_keys = sorted(f.keys(), key=lambda x: int(re.search(r'\d+', x).group()))
        selected_keys = step_keys[::stride]
        
        # extract all star mass and sound speed data
        for step_key in selected_keys:
            densities.append(np.array(f[step_key]['rho']))
            pressures.append(np.array(f[step_key]['p']))
            x_pos.append(np.array(f[step_key]['x']))
            y_pos.append(np.array(f[step_key]['y']))
            z_pos.append(np.array(f[step_key]['z']))
            times.append(np.array(f[step_key].attrs['time']))
            star_x.append(np.array(f[step_key].attrs['star::x']))
            star_y.append(np.array(f[step_key].attrs['star::y']))
            star_z.append(np.array(f[step_key].attrs['star::z']))

    return selected_keys, densities, pressures, x_pos, y_pos, z_pos, times, star_x, star_y, star_z

def calc_radius(x, y, z, star_x, star_y, star_z, ts, particle): 
    x_diff = x[ts][particle] - star_x[ts]
    y_diff = y[ts][particle] - star_y[ts]
    z_diff = z[ts][particle] - star_z[ts]

    r2 = (x_diff * x_diff) + (y_diff * y_diff)  + (z_diff * z_diff)
    return r2
    
def select_particles(timesteps, x, y, z, star_x, star_y, star_z): 
    particle_indices = np.array(timesteps)
    for step in range(0, len(timesteps), 1): 
        nr_particles = len(x[step])
        for p in range(0, nr_particles, 1):
            if calc_radius(x, y, z, star_x, star_y, star_z, step, p) <= threshold_radius: 
                particle_indices[step].append(p)



if __name__ == "__main__":
    ts, d, p, x, y, z, times, sx, sy, sz = read_hdf5_data(f5,1)
