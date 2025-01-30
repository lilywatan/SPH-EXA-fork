import vertical_profile as vp
import matplotlib.pyplot as plt
import h5py
import re
import seaborn as sns
from math import pi 
from scipy.spatial import cKDTree
from matplotlib.colors import Normalize
import numpy as np
plots = './output/plots/vertical-profile'

def plot_kin_visc(timesteps, c, r, z, m, disk_mask, stride, num_bins, alpha, beta):
    radial_bins = None
    plt.figure()
    for t in timesteps: 
        index = int(t/(10*stride))
        c_disk = c[index][disk_mask[index]]
        r_disk = r[index][disk_mask[index]]
        z_disk = z[index][disk_mask[index]]
        m_disk = m[index][disk_mask[index]]

        # Define radial bins only once
        if radial_bins is None:
            r_min, r_max = np.min(r_disk), np.max(r_disk)
            radial_bins = np.linspace(r_min, r_max, num_bins + 1)
            bin_centers = 0.5 * (radial_bins[:-1] + radial_bins[1:])

        # Compute RMS scale height for each radial bin
        rms_values = []
        nu_values = []
        for i in range(len(radial_bins) - 1):
            mask = (r_disk >= radial_bins[i]) & (r_disk < radial_bins[i+1])
            m_selected = m_disk[mask]
            z_selected = z_disk[mask]
            c_selected = c_disk[mask]
            
            rms = vp.scale_height_rms(m_selected, z_selected)
            rms_values.append(rms/bin_centers[i])
            nu_values.append(alpha * np.mean(c_selected) * rms)
        print(nu_values)
        plt.plot(bin_centers, nu_values, marker='o', label=f'Timestep {t}')
    plt.xlabel('Radius')
    plt.ylabel('kinematic viscosity')
    plt.title(f'kinematic viscosity of disk, momentum criterion, beta={beta}')
    plt.grid()
    plt.legend(loc='upper right', fontsize='small')
    fname = f'/kin-visc-rms/kin_visc_rms_100_mom_{beta}.png'
    plt.savefig(plots + fname)
    plt.show()
    

if __name__ == '__main__':
    stride=1
    ts, d, p, m, x, y, z, h, c_s, times, sx, sy, sz, sm = vp.read_hdf5_data(vp.run_100_beta,stride)
    r, disk_particles = vp.particle_radii2(x, y, z, m, sx, sy, sz, sm, ts)
    plot_kin_visc([1000, 50000, 10000], c_s, r, z, m, disk_particles, stride, 200, 0.232, "2pi")