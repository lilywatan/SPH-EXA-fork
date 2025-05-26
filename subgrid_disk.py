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

run_subgrid_star0 = '/home/lwatan/scratch/run_subgrid_beta_planet_star0_h5.hdf5'
run_subgrid_star0_lim12 = '/home/lwatan/scratch/run_subgrid_beta_planet_star0_1-2lim.hdf5'
run_subgrid_star0_lim34 = '/home/lwatan/scratch/run_subgrid_beta_planet_star0_3-4lim.hdf5'
run_subgrid_star0_radlim = '/home/lwatan/scratch/run_subgrid_beta_planet_star0_radial_lim_2.hdf5'
surface_density_min = None
surface_density_max = None
plots = '/home/lwatan/data/SPH-EXA-fork/output/plots/subgrid/radlim/'

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
    disk_r = []
    disk_m = []
    disk_sigma0 = []

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
            disk_r.append(np.array(f[step_key].attrs['disk::r']))
            disk_m.append(np.array(f[step_key].attrs['disk::m']))
            disk_sigma0.append(np.array(f[step_key].attrs['disk::sigma0']))



    return selected_keys, densities, pressures, masses, x_pos, y_pos, z_pos, h, c, times, star_x, star_y, star_z, star_m, disk_r0, disk_r, disk_m, disk_sigma0

def print_radii(r0):
    for i in r0:
        print(i, '\n')

def plot_subgrid_disk(t, x, y, disk_mask, stride, dr0, dr, star_x, star_y):
    print("length x_disk: ", len(x), "\n")
    index = int(t / (10 * stride))
    print("index: ", index, "\n")
    x_disk = x[index][disk_mask[index]]
    y_disk = y[index][disk_mask[index]]

    star_center = (star_x[index], star_y[index])
    circle_r0 = Circle(star_center, dr0[index], color='red', fill=False, linewidth=2, alpha=0.5)
    circle_r = Circle(star_center, dr[index], color='green', fill=False, linewidth=2)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.scatter(x_disk, y_disk,marker='.')
    ax.add_patch(circle_r)
    ax.add_patch(circle_r0)

    plt.title(f'Subgrid Disk Boundaries, (beta=2π)')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')

    fname = f'subgrid_star0_{t}_2π_corr.pdf'
    plt.savefig(plots + fname, bbox_inches='tight', dpi=300)
    plt.show()

def plot_radial_evolution(timesteps, dr0, dr):
    plt.figure(figsize=(8,6))
    for t in timesteps: 
        index = int(t / (1000 * stride))
        plt.plot(t, dr0[index])
        plt.plot(t, dr[index])

    plt.xlabel(r'Timestep')
    plt.ylabel("Radius [AU]")
    plt.title(f'Evolution of Subgrid Disk Radii over Time, beta=2π')
    plt.grid()
    plt.legend(loc='upper right', fontsize='small')
    fname = f'radial_evolution.pdf'
    plt.savefig(plots + fname)
    plt.show()

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
    print("length x_disk: ", len(x), "\n")
    index = int(t / (1000 * stride))
    x_disk = x[index][disk_mask[index]]
    y_disk = y[index][disk_mask[index]]
    
    # subgrid disk boundary
    star_center = (star_x, star_y)
    circle = Circle(star_center, r0, color='red', fill=False, linewidth=2)
    circle_r = Circle(star_center, r, color='red', fill=False, linewidth=2)

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))
    hb = ax.hexbin(x_disk, y_disk, C=surface_density, gridsize=grid_size, cmap='viridis', norm=norm)
    fig.colorbar(hb, ax=ax, label='Surface Density')
    ax.add_patch(circle)
    ax.add_patch(circle_r)
    
    # Title for fixed or varying scale
    scale_type = "Fixed" if fixed_scale else "Dynamic"
    plt.title(f'Surface Density at Timestep {t}, subgrid disk, (beta={beta}, {scale_type} Scale)')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Save the plot with scale type in filename
    fname = f'surface_density_1e6_subgrid_hexbin_{t}_{beta}_{scale_type.lower()}.png'
    plt.savefig(plots + fname, bbox_inches='tight', dpi=300)
    plt.show()

if __name__ == '__main__':
    stride=1
    ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm, d_r0, d_r, d_m, d_sig = read_hdf5_data_subgrid(run_subgrid_star0_radlim,stride)
    #print_radii(r0[800])
    r, disk_particles = surface_density.particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)

    min_step = 0
    step_interval = 100
    max_step = 1700

    for step in range(min_step, max_step, step_interval):
        # Compute surface density for the current step
        # sig = surface_density(step, d, x, y, h, m, disk_particles, stride)
        # surface_densities.append(sig)

        # # Compute global min/max for the fixed color scale
        # compute_global_min_max(surface_densities)

        # Generate plot for the current step with a fixed color scale
        #plot_surface_density(step, x, y, sig, disk_particles, stride, "2π", r0, sx, sy fixed_scale=True)
        plot_subgrid_disk(step, x, y, disk_particles, stride, d_r0, d_r, sx, sy)
