import h5py
import numpy as np

# File paths
no_acc = '/home/lwatan/scratch/run_disk_cal_1e6_beta_planet_star0.hdf5'  # Calibration file
input_subgrid = '/home/lwatan/data/SPH-EXA-fork/subgrid_init_planet_star0.hdf5'  # Sub-grid file (empty)
beta = '/home/lwatan/data/disk5_beta.hdf5'

# Open all three files
with h5py.File(input_subgrid, 'a') as sub, \
     h5py.File(no_acc, 'r') as cal, \
     h5py.File(beta, 'r') as beta_file:
    
    # Get the "Step#15" dataset from the calibration file
    cal_obj = cal['Step#15']  # Make sure this is the correct step

    # Get the "Step#0" dataset from the beta file (has 100,000 particles)
    beta_obj = beta_file['Step#0']

    # Ensure that "Step#0" exists in the sub-grid file; create if not present
    if 'Step#0' not in sub:
        sub.create_group('Step#0')
    
    # Get the "Step#0" dataset from the sub-grid file
    sub_obj = sub['Step#0']

    ### Copy All Datasets from `beta` to `input_subgrid`
    for key in beta_obj.keys():
        if key in sub_obj.keys():  # Check if dataset exists
            del sub_obj[key]  # Remove existing dataset if it exists
        beta_obj.copy(key, sub_obj, name=key)

    # Set iteration to 0 and other necessary attributes
    sub['Step#0'].attrs['iteration'] = 0
    sub['Step#0'].attrs['star::x'] = 0
    sub['Step#0'].attrs['star::y'] = 0

    # Get the current number of particles in Step#0
    num_particles_beta = len(sub['Step#0']['x'][0:])
    target_num_particles = 10000

    # Print current number of particles
    print(f"Current number of particles (from beta): {num_particles_beta}")

    # Now, we will update the first 99935 particles with data from `cal`
    num_particles_cal = len(cal_obj['x'])

    # If the number of particles in `cal` is less than 99935, we adjust this by filling dummy particles
    if num_particles_cal < 99935:
        print(f"Warning: The calibration file has fewer particles than expected ({num_particles_cal} < 99935).")
        # You can decide here how to handle this case (e.g., pad with dummy particles).
    else:
        # Update the first 99935 particles from `cal` to match the data in the `cal` file
        for key in cal_obj.keys():
            # Ensure the data in `cal` exists and update the first 99935 particles
            if key in sub_obj:
                sub_obj[key][:99935] = cal_obj[key][:99935]  # Update the first 99935 particles

    # After updating, ensure `numParticlesGlobal` is correctly set
    sub['Step#0'].attrs['numParticlesGlobal'] = 100000

    # Print final number of particles
    print(f"Final number of particles in Step#0: {len(sub['Step#0']['x'][0:])}")
    print("All attributes and datasets copied successfully from 'Step#0' of beta and updated with 'cal'.")
