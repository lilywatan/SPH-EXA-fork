import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import re
import seaborn as sns
from math import pi 
from scipy.spatial import cKDTree
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
import surface_density

run_subgrid_100 = '/home/lwatan/scratch/run_subgrid_100_beta.hdf5'
surface_density_min = None
surface_density_max = None

def read_hdf5_data_subgrid(file, stride=1):
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
    disk_r0 = []

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
            disk_r0.append(np.array(f[step_key].attrs['disk::r0']))



    return selected_keys, densities, pressures, masses, x_pos, y_pos, z_pos, h, c, times, star_x, star_y, star_z, star_m, disk_r0

def print_radii(r0):
    for i in r0:
        print(i, '\n')

def plot_surface_density(t, x, y, surface_density, disk_mask, stride, beta, r0, star_x, star_y, grid_size=200, fixed_scale=True): 
    # Global min/max for fixed scale
    global surface_density_min, surface_density_max

    # Determine normalization based on fixed or dynamic scale
    if fixed_scale:
        norm = Normalize(vmin=surface_density_min, vmax=surface_density_max)
    else:
        # Dynamic scale: Compute vmin and vmax for the current timestep
        current_min = surface_density.min()
        current_max = surface_density.max()
        norm = Normalize(vmin=current_min, vmax=current_max)

    # Index calculation for the given timestep
    index = int(t / (1 * stride))
    x_disk = x[index][disk_mask[index]]
    y_disk = y[index][disk_mask[index]]
    
    # subgrid disk boundary
    star_center = (star_x, star_y)
    circle = Circle(star_center, r0, color='red', fill=False, linewidth=2)

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))
    hb = ax.hexbin(x_disk, y_disk, C=surface_density, gridsize=grid_size, cmap='viridis', norm=norm)
    fig.colorbar(hb, ax=ax, label='Surface Density')
    ax.add_patch(circle)
    
    # Title for fixed or varying scale
    scale_type = "Fixed" if fixed_scale else "Dynamic"
    plt.title(f'Surface Density at Timestep {t}, subgrid disk, (beta={beta}, {scale_type} Scale)')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Save the plot with scale type in filename
    fname = f'surface_density_100_subgrid_hexbin_{t}_{beta}_{scale_type.lower()}.png'
    plt.savefig(surface_density.plots + fname, bbox_inches='tight', dpi=300)
    plt.show()

if __name__ == '__main__':
    stride=1
    ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm, r0 = read_hdf5_data_subgrid(run_subgrid_100,stride)
    print_radii(r0)
    r, disk_particles = surface_density.particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
    # Compute surface density for the specified timesteps
    sig_1 = surface_density.surface_density(1, d, x, y, h, m, disk_particles, stride)
    sig_10 = surface_density.surface_density(10, d, x, y, h, m, disk_particles, stride)
    sig_100 = surface_density.surface_density(100, d, x, y, h, m, disk_particles, stride)

    # Compute global min/max for the fixed color scale
    surface_density.compute_global_min_max([sig_100, sig_10, sig_1])

    # Generate plots with a fixed color scale
    plot_surface_density(1, x, y, sig_1, disk_particles, stride, r0, sx, sy, "2pi", fixed_scale=True)
    plot_surface_density(10, x, y, sig_10, disk_particles, stride, "2pi", fixed_scale=True)
    plot_surface_density(100, x, y, sig_100, disk_particles, stride, "2pi", fixed_scale=True)
