import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import re
import seaborn as sns
from math import pi 
from scipy.spatial import cKDTree
from matplotlib.colors import Normalize

run_20 = './output/run_disk_comb_20.hdf5'
run_20_beta = '/home/lwatan/data/SPH-EXA-fork/output/run_disk_comb_20_beta.hdf5'
plots = '/home/lwatan/data/SPH-EXA-fork/output/plots/'

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
            pressures.append(np.array(f[step_key]['p']))
            masses.append(np.array(f[step_key]['m']))
            x_pos.append(np.array(f[step_key]['x']))
            y_pos.append(np.array(f[step_key]['y']))
            z_pos.append(np.array(f[step_key]['z']))
            h.append(np.array(f[step_key]['h']))
            c.append(np.array(f[step_key]['c']))
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
                                 (y[step] - star_y[step])**2 + 
                                 (z[step] - star_z[step])**2)
        threshold_radius = 7.5
        step_disk_particles[:] = step_radii < threshold_radius

        radii.append(step_radii)
        disk_particles.append(step_disk_particles)
    return radii, disk_particles

# caluclate the scale height of the disk at radius r 
def calc_scale_height(c_s, r, m_s):
    angular_vel = np.sqrt(G * m_s / r**3)
    scale_height = c_s / angular_vel
    return scale_height

def plot_aspect_ratio(timesteps, c_s, r, m_s, disk_mask, stride, num_bins): 
    plt.figure()
    for t in timesteps: 
        index = int(t/(10*stride))
        c_disk = c_s[index][disk_mask[index]]
        r_disk = r[index][disk_mask[index]]
        m_s_disk = m_s[index]

        # bins for radii
        r_min, r_max = np.min(r_disk), np.max(r_disk)
        bins = np.linspace(r_min, r_max, num_bins + 1) 
        bin_centers = 0.5 * (bins[:-1] + bins[1:])  

        # average density in each bin
        average_ratio = np.zeros(num_bins)
        counts = np.zeros(num_bins) 

        for i in range(len(r_disk)):
            bin_index = np.digitize(r_disk[i], bins) - 1  # -1 to get the index in zero-based format
            if 0 <= bin_index < num_bins:  # check if bin_index is within valid range
                average_ratio[bin_index] += calc_scale_height(c_disk[i], r_disk[i], m_s_disk)/r_disk[i]
                counts[bin_index] += 1

        average_ratio /= counts
        average_ratio[counts == 0] = np.nan  # set bins with no particles to NaN for better plotting

        plt.plot(bin_centers, average_ratio, label=f'Timestep {t}', marker='.')
        plt.title('Aspect Ratio of Disk at Radius r')
        plt.xlabel('Radius')
        plt.ylabel('H(r)/r')
        plt.grid()
        plt.legend()
        fname = 'aspect_ratios.png'
        plt.savefig(plots + fname)  
        plt.show()

def plot_vertical_density(timesteps, z, densities, disk_mask, stride, num_bins=20):
    plt.figure()

    for t in timesteps:
        index = int(t/(stride*10))
        z_disk = z[index][disk_mask[index]]
        densities_disk = densities[index][disk_mask[index]]

        # binning 
        r_min, r_max = np.min(z_disk), np.max(z_disk)
        bins = np.linspace(r_min, r_max, num_bins + 1)  
        bin_centers = 0.5 * (bins[:-1] + bins[1:])  

        # average pressure per bin
        average_densities = np.zeros(num_bins)
        counts = np.zeros(num_bins)  

        for i in range(len(z_disk)):
            bin_index = np.digitize(z_disk[i], bins) - 1  # -1 to get the index in zero-based format
            if 0 <= bin_index < num_bins:  # check if bin_index is within valid range
                average_densities[bin_index] += densities_disk[i]
                counts[bin_index] += 1

        average_densities /= counts
        average_densities[counts == 0] = np.nan  # set bins with no particles to NaN for better plotting

        plt.plot(bin_centers,average_densities, label=f'Timestep {t}', marker='.')

        plt.title('Vertical Density Profile over Multiple Timesteps')
        plt.xlabel('Height')
        plt.ylabel('Density')
        plt.grid()
        plt.legend()
        fname = 'vertical_density.png'
        plt.savefig(plots + fname) 
        plt.show()

if __name__ == '__main__':
    stride=1
    ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm = read_hdf5_data(run_20,stride)
    r, disk_particles = particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
    plot_aspect_ratio([0, 5000, 10000, 15000, 20000], c_s, r, sm, disk_particles, stride, 100)
    plot_vertical_density([0, 5000, 10000, 15000, 20000], z, d, disk_particles, stride, 100)