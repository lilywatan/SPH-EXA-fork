import h5py

# File paths
no_acc = '/home/lwatan/data/SPH-EXA-fork/output/runs/run_disk_cal_1e6_beta_planet.hdf5'  # Calibration file
input_subgrid = '/home/lwatan/data/SPH-EXA-fork/subgrid_init_planet.hdf5'  # Sub-grid file (empty)
beta = '/home/lwatan/data/disk5_beta.hdf5'

# Open all three files
with h5py.File(input_subgrid, 'a') as sub, \
     h5py.File(no_acc, 'r') as cal, \
     h5py.File(beta, 'r') as beta_file:
    
    # Get the "Step#25" dataset from the calibration file
    cal_obj = cal['Step#16']
    
    # Get the "Step#0" dataset from the beta file
    beta_obj = beta_file['Step#0']

    # Ensure that "Step#0" exists in the sub-grid file; create if not present
    if 'Step#0' not in sub:
        sub.create_group('Step#0')
    
    # Get the "Step#0" dataset from the sub-grid file
    sub_obj = sub['Step#0']

    ### Copy only those Attributes that exist in the beta file.
    for attr in beta_obj.attrs.keys():
        if attr in cal_obj.attrs:
            sub_obj.attrs[attr] = cal_obj.attrs[attr]
        else:
            sub_obj.attrs[attr] = beta_obj.attrs[attr]

    ### Copy Datasets (Keys) from `no_acc` to `input_subgrid`
    for key in cal_obj.keys():
        if key in sub_obj:
            del sub_obj[key]  # Remove existing dataset if it exists to avoid conflicts
        cal_obj.copy(key, sub_obj, name=key)

    ### Copy Additional Datasets (Keys) from `beta` (if not present in `no_acc`)
    for key in beta_obj.keys():
        if key not in cal_obj:  # Only copy if it does not exist in `no_acc`
            if key in sub_obj:
                del sub_obj[key]  # Remove existing dataset if it exists
            beta_obj.copy(key, sub_obj, name=key)

    print("All attributes and datasets copied successfully from 'Step#25' to 'Step#0'.")
