import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import re
import seaborn as sns
from math import pi 
from scipy.spatial import cKDTree
from matplotlib.colors import Normalize

run_20 = './output/runs/run_disk_comb_20.hdf5'
run_20_beta = '/home/lwatan/data/SPH-EXA-fork/output/runs/run_disk_comb_20_beta.hdf5'
run_20_mom = './output/runs/run_disk_mom_20.hdf5'
run_20_mom_beta = './output/runs/run_disk_mom_20_beta.hdf5'
run_radial = './output/runs/run_disk_radial_20.hdf5'
run_50_beta = './output/runs/run_disk_mom_50_beta.hdf5'
run_5 = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/runs/run_disk_1J_5000.hdf5'
plots = './output/plots/vertical-profile'
run_100_beta = '/home/lwatan/scratch/run_disk_mom_100_beta.hdf5'
run_100_comb_beta = '/home/lwatan/scratch/run_disk_comb_100_beta.hdf5'
run_100_radial_beta = '/home/lwatan/scratch/run_disk_radial_100_beta.hdf5'
run_mom_1e6 = '/home/lwatan/scratch/run_disk_mom_1e6_beta.hdf5'
run_comb_1e6 = '/home/lwatan/scratch/run_disk_comb_1e6_2.hdf5'
run_rad_1e6 = '/home/lwatan/scratch/run_disk_radial_1e6_beta.hdf5'
run_comb_r1 = '/home/lwatan/scratch/run_disk_comb_1e6_r1_r1.hdf5'
run_half = '/home/lwatan/scratch/run_disk_half_1e6_25.hdf5'
run_double = '/home/lwatan/scratch/run_disk_double_1e6_short.hdf5'

# constants
G = 1.0

# extract data from hdf5 file with optional stride 
# returns arrays of selected keys (based on stride), density, pressure, x, y and z coords of particles & star
def read_hdf5_data(file, stride=1):
    densities = []
    pressures = []
    masses = []
    x_pos = []
    y_pos = []
    z_pos = []
    h = []
    c = []
    times = []
    star_x = []
    star_y = []
    star_z = []
    star_m = []

    with h5py.File(file, 'r') as f: 
        step_keys = sorted(f.keys(), key=lambda x: int(re.search(r'\d+', x).group()))
        selected_keys = step_keys[::stride]
        
        # extract all star mass and sound speed data
        for step_key in selected_keys:
            densities.append(np.array(f[step_key]['rho']))
            #pressures.append(np.array(f[step_key]['p']))
            masses.append(np.array(f[step_key]['m']))
            x_pos.append(np.array(f[step_key]['x']))
            y_pos.append(np.array(f[step_key]['y']))
            z_pos.append(np.array(f[step_key]['z']))
            h.append(np.array(f[step_key]['h']))
            #c.append(np.array(f[step_key]['c']))
            times.append(np.array(f[step_key].attrs['time']))
            star_x.append(np.array(f[step_key].attrs['star::x']))
            star_y.append(np.array(f[step_key].attrs['star::y']))
            star_z.append(np.array(f[step_key].attrs['star::z']))
            star_m.append(np.array(f[step_key].attrs['star::m']))

    return selected_keys, densities, pressures, masses, x_pos, y_pos, z_pos, h, c, times, star_x, star_y, star_z, star_m

# shortened version  
def particle_radii2(x, y, z, m, star_x, star_y, star_z, star_m, timesteps): 
    num_timesteps = len(timesteps)
    radii = []
    disk_particles = []

    for step in range(num_timesteps):
        nr_particles = len(x[step])  
        step_radii = np.zeros(nr_particles) 
        step_disk_particles = np.zeros(nr_particles, dtype=bool)
        
        # calculate the radii for each particle at this step
        step_radii[:] = np.sqrt((x[step] - star_x[step])**2 + 
                                 (y[step] - star_y[step])**2)
        threshold_radius = 50.0
        step_disk_particles[:] = step_radii < threshold_radius

        radii.append(step_radii)
        disk_particles.append(step_disk_particles)
    return radii, disk_particles

# caluclate the scale height of the disk at radius r 
def calc_scale_height(c_s, r, m_s):
    angular_vel = np.sqrt(G * m_s / r**3)
    scale_height = c_s / angular_vel
    return scale_height

def scale_height_rms(m, z, star_z): 
    if len(m) == 0:  # avoid division by zero for empty bins
        return 0
    dz = z - star_z; 
    return np.sqrt(np.sum(m * dz**2) / np.sum(m))

def get_scale_height_rms(timesteps, r, z, m, sz, disk_mask, stride, num_bins): 
    radial_bins = None
    scale_height_data = {}  # Store H(r) values per timestep
    aspect_ratio_data = {}  # Store H(r)/r values per timestep

    plt.figure()
    
    for t in timesteps: 
        index = int(t/(1000*stride))
        r_disk = r[index][disk_mask[index]]
        z_disk = z[index][disk_mask[index]]
        m_disk = m[index][disk_mask[index]]
        star_z = sz[index]

        # Define radial bins only once
        if radial_bins is None:
            r_min, r_max = np.min(r_disk), np.max(r_disk)
            radial_bins = np.linspace(r_min, r_max, num_bins + 1)
            bin_centers = 0.5 * (radial_bins[:-1] + radial_bins[1:])

        # Compute RMS scale height for each radial bin
        rms_values_H = []
        rms_values_ar = []
        for i in range(len(radial_bins) - 1):
            mask = (r_disk >= radial_bins[i]) & (r_disk < radial_bins[i+1])
            m_selected = m_disk[mask]
            z_selected = z_disk[mask]
            
            rms = scale_height_rms(m_selected, z_selected, star_z)
            rms_values_H.append(rms)
            rms_values_ar.append(rms / bin_centers[i] if bin_centers[i] != 0 else np.nan)  # Avoid division by zero
        
        # Store results for this timestep
        scale_height_data[t] = rms_values_H
        aspect_ratio_data[t] = rms_values_ar
    return bin_centers, scale_height_data, aspect_ratio_data 

def plot_all(t_double, t_rest, b_half, b_single, b_double, b_rad, s_half, s_single, s_double, s_rad, a_half, a_single, a_double, a_rad): 
    plt.figure(figsize=(8,6))
    plt.plot(b_double, s_double[0], marker='o', label=f'Initial {0}')
    for t in t_double: 
        plt.plot(b_double, s_double[t], marker='o', label=f'2.0 Momentum {t}')
    for t in t_rest: 
        plt.plot(b_half, s_half[t], marker='o', label=f'0.5 Momentum {t}')
        plt.plot(b_single, s_single[t], marker='o', label=f'1.0 Momentum {t}')
        plt.plot(b_rad, s_rad[t], marker='o', label=f'Radial Criterion {t}')
    plt.xlabel(r'Radius (r) $\left[ AU \right]$')
    plt.ylabel(r"Scale Height (H) $\left[ AU \right]$")
    plt.title(f'Scale Height vs Radius of Different Accretion Criteria')
    plt.grid()
    plt.legend(loc='upper right', fontsize='small')
    fname = f'/scale-height-rms/scale_height_rms_1e6_comparison.pdf'
    plt.savefig(plots + fname)
    plt.show()

    # aspect ratio
    plt.figure(figsize=(8,6))
    plt.plot(b_double, s_double[0], marker='o', label=f'Initial {0}')
    for t in t_double: 
        plt.plot(b_double, a_double[t], marker='o', label=f'2.0 Momentum {t}')
    for t in t_rest: 
        plt.plot(b_half, a_half[t], marker='o', label=f'0.5 Momentum {t}')
        plt.plot(b_single, a_single[t], marker='o', label=f'1.0 Momentum {t}')
        plt.plot(b_rad, a_rad[t], marker='o', label=f'Radial Criterion {t}')
    plt.xlabel(r'Radius (r) $\left[ AU \right]$')
    plt.ylabel("Aspect Ratio (H(r)/r)")
    plt.title(f'Aspect Ratio vs Radius of Different Accretion Criteria')
    plt.grid()
    plt.legend(loc='upper right', fontsize='small')
    fname = f'/ar-rms/ar_rms_1e6_comparison.pdf'
    plt.savefig(plots + fname)
    plt.show()

def gen_data(run_value, steps): 
    ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm = read_hdf5_data(run_value,stride)
    r, disk_particles = particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
    b, s, a = get_scale_height_rms(steps, r, z, m, sz, disk_particles, stride, 20)
    return b, s, a

if __name__ == '__main__':
    stride=1

    steps_rest = [7000, 100000, 250000]
    steps_double = [7000]

    b_half, s_half, a_half = gen_data(run_half, steps_rest)
    b_double, s_double, a_double = gen_data(run_double, steps_double)
    b_single, s_single, a_single = gen_data(run_mom_1e6, steps_rest)
    b_rad, s_rad, a_rad = gen_data(run_rad_1e6, steps_rest)

    plot_all(steps_double, steps_rest, b_half, b_single, b_double, b_rad, s_half, s_single, s_double, s_rad, a_half, a_single, a_double, a_rad)


