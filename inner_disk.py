import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import re

f5 = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_1J_5000.hdf5'
f5_double = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_2J_5000.hdf5'
f5_half = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_05J_5000.hdf5'
plots = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/plots/'

# constants
G = 1 
G2 = G * G

# extract data from hdf5 file with optional stride 
# returns arrays of selected keys (based on stride), density, pressure, x, y and z coords of particles & star
def read_hdf5_data(file, stride=1):
    densities = []
    pressures = []
    masses = []
    x_pos = []
    y_pos = []
    z_pos = []
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
            times.append(np.array(f[step_key].attrs['time']))
            star_x.append(np.array(f[step_key].attrs['star::x']))
            star_y.append(np.array(f[step_key].attrs['star::y']))
            star_z.append(np.array(f[step_key].attrs['star::z']))
            star_m.append(np.array(f[step_key].attrs['star::m']))

    return selected_keys, densities, pressures, masses, x_pos, y_pos, z_pos, times, star_x, star_y, star_z, star_m

# return the squared radial distance between particle & star
def calc_radius(x, y, z, star_x, star_y, star_z, ts, particle): 
    x_diff = x[ts][particle] - star_x[ts]
    y_diff = y[ts][particle] - star_y[ts]
    z_diff = z[ts][particle] - star_z[ts]

    r2 = (x_diff * x_diff) + (y_diff * y_diff)  + (z_diff * z_diff)
    return r2

# calculates radial distance from star for each particle at each timestep 
# returns array of radii 
def particle_radii(x, y, z, star_x, star_y, star_z, timesteps): 
    radii = np.array(timesteps)
    for step in range(0, len(timesteps), 1):
        nr_particles = len(x[step])
        for p in range(0, nr_particles, 1): 
            r2 = calc_radius(x, y, z, star_x, star_y, star_z, step, p)
            radii[step][p] = (np.sqrt(r2))
    return radii 

# shortened version  
def particle_radii2(x, y, z, m, star_x, star_y, star_z, star_m, timesteps): 
    num_timesteps = len(timesteps)
    # Initialize lists to hold the arrays for each time step
    radii = []
    disk_particles = []

    for step in range(num_timesteps):
        nr_particles = len(x[step])  # Number of particles at this time step
        # Create numpy arrays for this time step
        step_radii = np.zeros(nr_particles)  # Initialize with zeros
        step_disk_particles = np.zeros(nr_particles, dtype=bool)  # Boolean array for selection
        
        # Calculate the radii for each particle at this step
        step_radii[:] = np.sqrt((x[step] - star_x[step])**2 + 
                                 (y[step] - star_y[step])**2 + 
                                 (z[step] - star_z[step])**2)

        # Mark disk particles based on the threshold
        threshold_radius = step_radii * np.cbrt(m[step]/(3 * star_m[step]))
        step_disk_particles[:] = step_radii < threshold_radius

        # Append the results for this step to the lists
        radii.append(step_radii)
        disk_particles.append(step_disk_particles)
    print("len of radii: ", len(radii), " len of disk_particles: ", len(disk_particles))
    return radii, disk_particles

# iterate through timesteps & particle arrays to select the particles in inner disk 
# threshold radius set based on 
def select_particles(timesteps, x, y, z, star_x, star_y, star_z, r): 
    disk_particles = np.array(timesteps)
    for step in range(0, len(timesteps), 1): 
        nr_particles = len(x[step])
        for p in range(0, nr_particles, 1):
            if r[step][p] <= threshold_radius: 
                disk_particles[step].append(p)
    return disk_particles

# def surface_density(timesteps, density, r, z): 

def plot_density(t, x, y, z, densities, disk_mask):
    x_disk = x[t][disk_mask[t]]
    y_disk = y[t][disk_mask[t]]
    z_disk = z[t][disk_mask[t]]
    densities_disk = densities[t][disk_mask[t]]
    print("Filtered X values:", x_disk)
    print("Filtered Y values:", y_disk)
    print("Filtered Z values:", z_disk)
    print("Filtered Densities:", densities_disk)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x_disk, y_disk, z_disk, c=densities_disk, cmap='viridis')
    plt.colorbar(scatter, label='Density')
    ax.set_title(f'Density Distribution in 3D at timestep {t}')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_xlim([np.min(x_disk), np.max(x_disk)])
    ax.set_ylim([np.min(y_disk), np.max(y_disk)])
    ax.set_zlim([np.min(z_disk), np.max(z_disk)])
    fname = f'density_distribution_1J.png'
    plt.savefig(plots + fname)
    plt.show()

def plot_pressure(t, x, y, z, pressures, disk_mask):
    x_disk = x[t][disk_mask[t]]
    y_disk = y[t][disk_mask[t]]
    z_disk = z[t][disk_mask[t]]
    pressures_disk = pressures[t][disk_mask[t]]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x_disk, y_disk, z_disk, c=pressures_disk, cmap='viridis')
    plt.colorbar(scatter, label='Pressure')
    ax.set_title(f'Pressure Distribution in 3D at timestep {t}')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_zlabel('Z Position')
    ax.set_xlim([np.min(x_disk), np.max(x_disk)])
    ax.set_ylim([np.min(y_disk), np.max(y_disk)])
    ax.set_zlim([np.min(z_disk), np.max(z_disk)])
    fname = f'pressure_distribution_1J.png'
    plt.savefig(plots + fname)
    plt.show()


if __name__ == "__main__":
    ts, d, p, m, x, y, z, times, sx, sy, sz, sm = read_hdf5_data(f5,1)
    r, disk_particles = particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
    plot_density(1, x, y, z, d, disk_particles)
    plot_pressure(1, x, y, z, p, disk_particles)

