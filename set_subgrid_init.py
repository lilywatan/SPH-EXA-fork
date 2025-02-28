import h5py

# File paths
no_acc = '/home/lwatan/scratch/run_disk_cal_1e6_beta_planet_star0.hdf5'  # Calibration file
input_subgrid = '/home/lwatan/data/SPH-EXA-fork/subgrid_init_planet_star0.hdf5'  # Sub-grid file (empty)
beta = '/home/lwatan/data/disk5_beta.hdf5'

# Open all three files
with h5py.File(input_subgrid, 'a') as sub, \
     h5py.File(no_acc, 'r') as cal, \
     h5py.File(beta, 'r') as beta_file:
    
    # Get the "Step#25" dataset from the calibration file (was mistakenly Step#16)
    cal_obj = cal['Step#16']  # Make sure this is the correct step

    # Get the "Step#0" dataset from the beta file
    beta_obj = beta_file['Step#0']

    # Ensure that "Step#0" exists in the sub-grid file; create if not present
    if 'Step#0' not in sub:
        sub.create_group('Step#0')
    
    # Get the "Step#0" dataset from the sub-grid file
    sub_obj = sub['Step#0']

    ### Copy only those Attributes that exist in the beta file.
    for attr in cal_obj.attrs.keys():
        if attr in sub_obj.attrs:  # Fixed: Check in attrs, not sub_obj directly
            del sub_obj.attrs[attr]  # Remove existing dataset if it exists
        sub_obj.attrs[attr] = cal_obj.attrs[attr]

    ### Copy Additional attributes from `beta` (if not present in `no_acc`)
    for attr in beta_obj.attrs.keys():
        if attr not in cal_obj.attrs:  # Only copy if it does not exist in `no_acc`
            if attr in sub_obj.attrs:  # Fixed: Check in attrs
                del sub_obj.attrs[attr]  # Remove existing dataset if it exists
            sub_obj.attrs[attr] = beta_obj.attrs[attr]

    ### Copy Datasets (Keys) from `no_acc` to `input_subgrid`
    for key in cal_obj.keys():
        if key in sub_obj.keys():  # Fixed: Check in keys()
            del sub_obj[key]  # Remove existing dataset if it exists
        cal_obj.copy(key, sub_obj, name=key)

    ### Copy Additional Datasets (Keys) from `beta` (if not present in `no_acc`)
    for key in beta_obj.keys():
        if key not in cal_obj:  # Only copy if it does not exist in `no_acc`
            if key in sub_obj.keys():  # Fixed: Check in keys()
                del sub_obj[key]  # Remove existing dataset if it exists
            beta_obj.copy(key, sub_obj, name=key)

    # Set iteration to 0
    sub['Step#0'].attrs['iteration'] = 0
    sub['Step#0'].attrs['star::x'] = 0
    sub['Step#0'].attrs['star::y'] = 0
    print(sub['Step#0']['x'].shape)
    # sub['Step#0'].attrs['star::removal_limit_h'] 
    del sub['Step#0'].attrs['star::inner_size'] 


    print("All attributes and datasets copied successfully from 'Step#25' to 'Step#0'.")
