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


run_20 = './output/runs/run_disk_rad_20.hdf5'
run_20_beta = '/home/lwatan/data/SPH-EXA-fork/output/runs/run_disk_rad_20_beta.hdf5'
run_J_20 = './output/runs/run_disk_mom_20.hdf5'
run_J_20_beta = './output/runs/run_disk_mom_20_beta.hdf5'
run_radial = './output/runs/run_disk_radial_20.hdf5'
plots = '/home/lwatan/data/SPH-EXA-fork/output/plots/'
run_50_beta = './output/runs/run_disk_mom_50_beta.hdf5'
run_100_beta = '/home/lwatan/scratch/run_disk_mom_100_beta.hdf5'
run_100_radial_beta = '/home/lwatan/scratch/run_disk_radial_100_beta.hdf5'
run_500_beta = '/home/lwatan/scratch/run_disk_mom_500_beta.hdf5'
run_1e6_beta = '/home/lwatan/scratch/run_disk_mom_1e6_beta.hdf5'
run_calibration = './output/runs/run_disk_cal_1e6_beta_2.hdf5'
run_cal_planet = '/home/lwatan/scratch/run_disk_cal_1e6_beta_planet_3.hdf5'

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
            #densities.append(np.array(f[step_key]['rho']))
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
                                 (y[step] - star_y[step])**2 )
        threshold_radius = 7.5
        step_disk_particles[:] = step_radii < threshold_radius

        radii.append(step_radii)
        disk_particles.append(step_disk_particles)
    return radii, disk_particles

# Kernel function according to paper
def W(r, h, sigma=10/(7*pi)):
    q = r/h 
    if 0 <= q <= 1: 
        return (sigma / h**2) * (1 - 1.5 * q**2 + 0.75 * q**3)
    elif 1 < q <= 2: 
        return (sigma / h**2) * 0.25*(2 - q)**3
    else: 
        return 0

# calculate surface density at timestep t for all particles 
def surface_density(t, rho, x, y, h, m, disk_mask, stride): 
    index = int(t/(1000*stride))
    print(len(rho))
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
    
def plot_surface_density(t, x, y, surface_density, disk_mask, stride, beta, grid_size=200, fixed_scale=True): 
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
    index = int(t / (1000 * stride))
    x_disk = x[index][disk_mask[index]]
    y_disk = y[index][disk_mask[index]]
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.hexbin(x_disk, y_disk, C=surface_density, gridsize=grid_size, cmap='viridis', norm=norm)
    plt.colorbar(label='Surface Density')
    
    # Title for fixed or varying scale
    scale_type = "Fixed" if fixed_scale else "Dynamic"
    plt.title(f'Surface Density at Timestep {t}, no accretion (beta={beta}, {scale_type} Scale)')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Save the plot with scale type in filename
    fname = f'surface_density_1e6_cal_hexbin_{t}_{beta}_{scale_type.lower()}.png'
    plt.savefig(plots + fname, bbox_inches='tight', dpi=300)
    plt.show()

def plot_particles(t, x, y, h, disk_mask, stride, beta="2pi"):
    index = int(t / (1000 * stride))
    
    # Apply disk mask to extract particle positions and smoothing lengths
    x_disk = x[index][disk_mask[index]]
    y_disk = y[index][disk_mask[index]]
    h_disk = h[index][disk_mask[index]]
    
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    
    # Plot the particle positions
    plt.scatter(x_disk, y_disk, marker='.', label='Particles')
    
    # Plot circles for each particle using its own smoothing length
    # Only add labels once for the legend (for the first instance of each)
    #first_2h = True
    #first_3h = True
    #for xi, yi, hi in zip(x_disk, y_disk, h_disk):
    #    circle_2h = Circle(
    #        (xi, yi), 2 * hi, 
    #        color='red', fill=False, linestyle='--', linewidth=1,
    #        label='2h' if first_2h else None
    #    )
    #    circle_3h = Circle(
    #        (xi, yi), 3 * hi, 
    #        color='blue', fill=False, linestyle='--', linewidth=1,
    #        label='3h' if first_3h else None
    #    )
    #    ax.add_patch(circle_2h)
    #    ax.add_patch(circle_3h)
    #    first_2h = False
    #   first_3h = False
    
    # Set title and axis labels
    plt.title(f'Particles at timestep {t}, no accretion (beta={beta})')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Add legend
    plt.legend()
    
    # Save the plot (make sure 'plots' variable is defined with a valid directory path)
    fname = f'particles_1e6_cal_planet_3_{t}_{beta}.png'
    plt.savefig(plots + fname, bbox_inches='tight', dpi=300)
    plt.show()

    
if __name__ == '__main__':
    stride=1
    ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm = read_hdf5_data(run_cal_planet,stride)
    r, disk_particles = particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
    # Compute surface density for the specified timesteps
    # Define the step interval
    step_interval = 10000
    max_step = 70000  # Adjust this to the maximum step in your simulation

    # Prepare an empty list to store surface density arrays for computing global min/max
    surface_densities = []

    # Loop through steps in increments of 100,000
    #for step in range(step_interval, max_step + step_interval, step_interval):
        # Compute surface density for the current step
        #sig = surface_density(step, d, x, y, h, m, disk_particles, stride)
        #surface_densities.append(sig)

    # Compute global min/max for the fixed color scale
    #compute_global_min_max(surface_densities)

    # Loop again to generate plots after computing global min/max
    for step in range(10000, max_step + step_interval, step_interval):
        # Generate plot for the current step with a fixed color scale
        #plot_surface_density(step, x, y, sig, disk_particles, stride, "2pi", fixed_scale=True)
        plot_particles(step, x, y, h, disk_particles, stride, "2pi")


    # Generate plots with a dynamic color scale
    #plot_surface_density(100000, x, y, sig_100k, disk_particles, stride, "2pi", fixed_scale=False)
    #plot_surface_density(75000, x, y, sig_75k, disk_particles, stride, "2pi", fixed_scale=False)
    #plot_surface_density(50000, x, y, sig_50k, disk_particles, stride, "2pi", fixed_scale=False)
    #plot_surface_density(25000, x, y, sig_25k, disk_particles, stride, "2pi", fixed_scale=False)
    #plot_surface_density(10, x, y, sig_10, disk_particles, stride, "2pi", fixed_scale=False)


