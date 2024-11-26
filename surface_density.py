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

# Kernel function according to paper
def W(r, h, sigma=10/7*pi):
    q = r/h 
    if 0 <= q <= 1: 
        return (sigma / h**2) * (1 - 1.5 * q**2 + 0.75 * q**3)
    elif 1 < q <= 2: 
        return (sigma / h**2) * 0.25*(2 - q)**3
    else: 
        return 0

# calculate surface density at timestep t for all particles 
def surface_density(t, rho, x, y, h, m, disk_mask, stride): 
    index = int(t/(10*stride))
    rho_disk = rho[index][disk_mask[index]]
    h_disk = h[index][disk_mask[index]]
    x_disk = x[index][disk_mask[index]]
    y_disk = y[index][disk_mask[index]]
    m_disk = m[index][disk_mask[index]]
    num_particles = len(x_disk)

    surface_density = np.zeros(num_particles)
    positions = np.column_stack((x_disk, y_disk))
    tree = cKDTree(positions)

    for j in range(1, num_particles):
        density = 0.0

        indices = tree.query_ball_point(positions[j], h_disk[j])
        for i in indices: 
            if i != j: 
                r_ij = np.linalg.norm(positions[j] - positions[i])
                W_ij = W(r_ij, h_disk[i])
                density += m_disk[i] * W_ij
        surface_density[j] = density
    return surface_density

surface_density_min = None
surface_density_max = None

def compute_global_min_max(surface_densities):
    global global_surface_density_min, global_surface_density_max
    # Find the global min/max across all surface densities
    all_values = np.concatenate([sd.flatten() for sd in surface_densities])
    global_surface_density_min = np.min(all_values)
    global_surface_density_max = np.max(all_values)
    
def plot_surface_density(t, x, y, surface_density, disk_mask, stride, grid_size=200): 
    global surface_density_min, surface_density_max
    
    norm = Normalize(vmin=global_surface_density_min, vmax=global_surface_density_max)
    index = int(t/(10*stride))

    x_disk = x[index][disk_mask[index]]
    y_disk = y[index][disk_mask[index]]
    # Scatter plot, with color representing the surface density
    plt.figure(figsize=(8, 6))
    #scatter = plt.scatter(x_disk, y_disk, c=surface_density, cmap="viridis", s=20, norm=norm)
    plt.hexbin(x_disk, y_disk, C=surface_density, gridsize=200, cmap='viridis', norm=norm)
    plt.hist2d(x_disk, y_disk, weights=surface_density, bins=grid_size, cmap='viridis', norm=norm)
    plt.colorbar(label='Surface Density')  # Show color scale
    plt.title(f'Surface Density at Timestep {t} (beta=2pi)')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    fname = f'surface density_20_beta_{t}.png'
    plt.savefig(plots + fname)
    plt.show()
    
if __name__ == '__main__':
    stride=1
    ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm = read_hdf5_data(run_20_beta,stride)
    r, disk_particles = particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
    sig_20k = surface_density(20000, d, x, y, h, m, disk_particles, stride)
    sig_15k = surface_density(15000, d, x, y, h, m, disk_particles, stride)
    sig_10k = surface_density(10000, d, x, y, h, m, disk_particles, stride)
    sig_5k = surface_density(5000, d, x, y, h, m, disk_particles, stride)
    sig_10 = surface_density(10, d, x, y, h, m, disk_particles, stride)

    compute_global_min_max([sig_20k, sig_15k, sig_10k, sig_5k, sig_10])

    plot_surface_density(20000, x, y, sig_20k, disk_particles, stride)
    plot_surface_density(15000, x, y, sig_15k, disk_particles, stride)
    plot_surface_density(10000, x, y, sig_10k, disk_particles, stride)
    plot_surface_density(5000, x, y, sig_5k, disk_particles, stride)
    plot_surface_density(10, x, y, sig_10, disk_particles, stride)