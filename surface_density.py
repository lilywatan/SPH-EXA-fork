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
plots = '/home/lwatan/data/SPH-EXA-fork/output/plots/surface-density/'
run_50_beta = './output/runs/run_disk_mom_50_beta.hdf5'
run_100_beta = '/home/lwatan/scratch/run_disk_mom_100_beta.hdf5'
run_100_radial_beta = '/home/lwatan/scratch/run_disk_radial_100_beta.hdf5'
run_500_beta = '/home/lwatan/scratch/run_disk_mom_500_beta.hdf5'
run_1e6_beta = '/home/lwatan/scratch/run_disk_comb_r1_beta.hdf5'
run_1e6_beta = '/home/lwatan/scratch/run_disk_comb_r1_beta.hdf5'
run_calibration = './output/runs/run_disk_cal_1e6_beta_2.hdf5'
run_cal_planet_3 = '/home/lwatan/scratch/run_disk_cal_1e6_beta_planet_3.hdf5'
run_cal_planet_2 = '/home/lwatan/scratch/run_disk_cal_1e6_beta_planet_2.hdf5'
run_comb_r1 = '/home/lwatan/scratch/run_disk_comb_r1_beta.hdf5'
#run_half = '/home/lwatan/scratch/run_disk_half_2.hdf5'
run_double = '/home/lwatan/scratch/run_disk_radial_1e6_beta.hdf5'
# calibration runs
run_cal_star0 = '/home/lwatan/scratch/run_disk_cal_1e6_beta_planet_star0.hdf5'
run_cal_no_hlim = '/home/lwatan/scratch/run_disk_cal_1e6_beta_no_hlim.hdf5'
run_cal_h5 = '/home/lwatan/scratch/run_disk_cal_1e6_beta_h5.hdf5'

# more runs
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
        threshold_radius = 50.0
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
    print(len(m), " ", index)
    #rho_disk = rho[index][disk_mask[index]]
    h_disk = h[index][disk_mask[index]]
    x_disk = x[index][disk_mask[index]]
    y_disk = y[index][disk_mask[index]]
    m_disk = m[index][disk_mask[index]]
    num_particles = len(x_disk)

    surface_density = np.zeros(num_particles)
    positions = np.column_stack((x_disk, y_disk))
    tree = cKDTree(positions)

    print("num_particles: ", num_particles)
    for j in range(num_particles):
        density = 0.0

        indices = tree.query_ball_point(positions[j], h_disk[j])
        for i in indices: 
            if i != j: 
                r_ij = np.linalg.norm(positions[j] - positions[i])
                W_ij = W(r_ij, h_disk[i])
                density += m_disk[i] * W_ij
        surface_density[j] = density
    print("length sd: ", len(surface_density))
    return surface_density

surface_density_min = None
surface_density_max = None

def compute_global_min_max(surface_densities):
    global global_surface_density_min, global_surface_density_max
    # Find the global min/max across all surface densities
    all_values = np.concatenate([sd.flatten() for sd in surface_densities])
    global_surface_density_min = np.min(all_values)
    global_surface_density_max = np.max(all_values)
    
def plot_surface_density(t, x, y, surface_density, disk_mask, stride, beta, run, grid_size=200, fixed_scale=True): 
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
    print(x_disk.shape, " ", y_disk.shape, " ", surface_density.shape)
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.hexbin(x_disk, y_disk, C=surface_density, gridsize=grid_size, cmap='viridis', norm=norm)
    plt.colorbar(label='Surface Density')
    
    # Title for fixed or varying scale
    scale_type = "Fixed" if fixed_scale else "Dynamic"
    plt.title(rf"Surface Density at Timestep {t} $\frac{{\left[ M_{{\odot}} \right]}}{{\left[ AU^2 \right]}}$, (β={beta})")
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Save the plot with scale type in filename
    fname = f'surface_density_1e6_{run}_hexbin_{t}_{beta}_{scale_type.lower()}_25.pdf'
    plt.savefig(plots + fname, bbox_inches='tight', dpi=300)
    plt.show()

def plot_particles(t, x, y, h, disk_mask, stride, beta, run):
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
    plt.title(f'Particles at Timestep {t}, No Accretion (β={beta})')
    plt.xlabel('X Position')
    plt.ylabel('Y Position')
    
    # Add legend
    plt.legend()
    
    # Save the plot (make sure 'plots' variable is defined with a valid directory path)
    fname = f'particles_1e6_{run}_{t}_{beta}_r1.pdf'
    plt.savefig(plots + fname, bbox_inches='tight', dpi=300)
    plt.show()

    
if __name__ == '__main__':
    stride=1
    # runs = {
    # #"run_comb_r1": run_comb_r1,
    # "run_half": run_half,
    # #"run_double": run_double
    # }   
    # for run_name, run_value in runs.items(): 
    #     print(f"starting with run {run_name}")
    #     ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm = read_hdf5_data(run_value,stride)
    #     r, disk_particles = particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
    #     run_type = run_name.split('_')[1]
    #     # Compute surface density for the specified timesteps
    #     # Define the step interval
    #     min_step = 0
    #     step_interval = 50000
    #     max_step = 260000  # Adjust this to the maximum step in your simulation

    #     # Prepare an empty list to store surface density arrays for computing global min/max
    #     surface_densities = []

    #     # Loop through steps in increments of 100,000
    #     for step in range(min_step, max_step, step_interval):
    #         # Compute surface density for the current step
    #         sig = surface_density(step, d, x, y, h, m, disk_particles, stride)
    #         surface_densities.append(sig)

    #         # Compute global min/max for the fixed color scale
    #         compute_global_min_max(surface_densities)

    #         # Generate plot for the current step with a fixed color scale
    #         plot_surface_density(step, x, y, sig, disk_particles, stride, "2π", run_type, fixed_scale=True)

    #     print(f"finished with run {run_name}")
    #     #plot_particles(step, x, y, h, disk_particles, stride, "2pi")

    cal = {
    "run_cal_star0": run_cal_star0,
    #"run_cal_no_hlim": run_cal_no_hlim,
    #"run_cal_h5": run_cal_h5
    }   

    for run_name, run_value in cal.items(): 
        print(f"starting with run {run_name}")
        ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm = read_hdf5_data(run_value,stride)
        r, disk_particles = particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
        run_type = run_name.split('_')[1]

        min_step = 0
        step_interval = 50000
        max_step = 260000

        for step in range(min_step, max_step, step_interval):
            plot_particles(step, x, y, h, disk_particles, stride, "2π", run_type)

        print(f"finished with run {run_name}")


    # Generate plots with a dynamic color scale
    #plot_surface_density(100000, x, y, sig_100k, disk_particles, stride, "2pi", fixed_scale=False)
    #plot_surface_density(75000, x, y, sig_75k, disk_particles, stride, "2pi", fixed_scale=False)
    #plot_surface_density(50000, x, y, sig_50k, disk_particles, stride, "2pi", fixed_scale=False)
    #plot_surface_density(25000, x, y, sig_25k, disk_particles, stride, "2pi", fixed_scale=False)
    #plot_surface_density(10, x, y, sig_10, disk_particles, stride, "2pi", fixed_scale=False)


