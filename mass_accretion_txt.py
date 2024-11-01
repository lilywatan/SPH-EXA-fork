# script to plot mass accretion rate from simulation vs theoretical mass accretion rate based on bondi accretion
# (G / c_s^3)

import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
from scipy.interpolate import CubicSpline 
from scipy.interpolate import UnivariateSpline as us
from scipy.stats import zscore
import re

file = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/star_mass_sound_speed_1000.txt'
f5 = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_sound_speed_1000.hdf5'
f5_double = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_2J.hdf5'
f5_half = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_05J.hdf5'
output = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/data_analysis_theoretical_1000.txt'
plots = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/plots/'

# constants
G = 1 
G2 = G * G
alpha = 4e-7 # scaling parameter for theoretical mass accretion rate

# extracts sound speed, timesteps & star mass from .txt file
# returns lists of masses, sound speeds & time steps
def read_txt_data(file, i=1):
    mass_list = []
    time_steps = []
    sound_speed = []
    with open(file, 'r') as f:
        # extract data from .txt file
        l = f.readlines()
        for i in range(0, len(l), i):
            # split the line by comma and strip whitespace
            c_s, time_step, mass_value = map(float, l[i].strip().split(','))
            sound_speed.append(c_s)
            time_steps.append(time_step)  # time difference per step (minDt)
            mass_list.append(mass_value)
    print(len(time_steps), "\n")
    return mass_list, sound_speed, time_steps

def read_hdf5_data(file, output=None, i=1):
    f5_mass_list = []
    f5_sound_speed = []
    with h5py.File(file, 'r') as f:
        step_keys = list(f.keys())
        for j in range(0, len(step_keys), i):
            star_mass_tmp = 0.0
            sound_speed_tmp = 0.0
            
            # take the average of 10 data points
            for k in range (0, i, 1):
                step_key = 'Step#' + str(i + k) 
                star_mass_tmp += f[step_key].attrs['star::m']
                sound_speed_tmp += sum(f[step_key]['c'])/len(f[step_key]['c'])
                if output != None: 
                    with open(output, 'a') as o: 
                        print(f[step_key].attrs['star::m'], ", ", file=o)
                        print(sum(f[step_key]['c'])/len(f[step_key]['c']), "\n", file=o)
            f5_mass_list.append(star_mass_tmp/10.0)
            f5_sound_speed.append(sound_speed_tmp/10.0)
        
    return f5_mass_list, f5_sound_speed

def read_hdf5_data_2(file, output=None, block_size=1):
    f5_mass_list = []
    f5_sound_speed = []

    with h5py.File(file, 'r') as f: 
        step_keys = sorted(f.keys(), key=lambda x: int(re.search(r'\d+', x).group()))

        # Initialize lists for temporary data storage
        star_mass_data = []
        sound_speed_data = []
        
        # Extract all star mass and sound speed data
        for step_key in step_keys:
            star_mass_data.append(f[step_key].attrs['star::m'])
            sound_speed_data.append(np.mean(f[step_key]['c']))

        # Convert lists to numpy arrays for vectorized operations
        star_mass_data = np.array(star_mass_data)
        sound_speed_data = np.array(sound_speed_data)

        # Determine number of full chunks
        num_chunks = len(star_mass_data) // block_size

        # Reshape arrays into chunks of 10 for averaging
        f5_mass_list = np.mean(star_mass_data[:num_chunks * block_size].reshape(-1, block_size), axis=1)
        f5_sound_speed = np.mean(sound_speed_data[:num_chunks * block_size].reshape(-1, block_size), axis=1)

        # Write to output if specified
        if output is not None:
            with open(output, 'a') as o:
                for mass, speed in zip(f5_mass_list, f5_sound_speed):
                    o.write(f"{mass}, {speed}\n")

        print(len(f5_mass_list), " ", len(f5_sound_speed))

    return f5_mass_list, f5_sound_speed


# calculate mass accretion rate from simulation data
# input: lists of masses & timesteps, return list of mass accretion rates
def calc_mass_accretion(masses, timesteps):
    mass_accretion_rate = []
    for i in range(1, len(masses)):
        mass_rate = (masses[i] - masses[i - 1]) / timesteps[i]
        # print(mass_rate, ",", time_steps[i], ",", mass_rate/time_steps[i], "\n", file=o)
        mass_accretion_rate.append(mass_rate)
    return mass_accretion_rate

# calculates the theoretical mass accretion rate based on bondi accretion
# input: list of sound speeds, return: list of theoretical accretion values
def calc_theoretical_accretion(c_s):
    theoretical_rate = []
    for i in range(1, len(c_s)):
        theoretical_rate.append(alpha * (G / (c_s[i] * c_s[i] * c_s[i]))) # M_dot = a * (G / c_s^3)
    return theoretical_rate

def filter_data(accretion_rates, theoretical_rates):
    accretion_rates = np.array(accretion_rates)
    theoretical_rates = np.array(theoretical_rates)

    # Compute Z-scores for accretion rates
    z_scores_accretion = zscore(accretion_rates)
    z_scores_theoretical = zscore(theoretical_rates)

    # Define a threshold for Z-score to detect outliers (e.g., 3 standard deviations)
    threshold = 2

    # Filter out outliers
    filtered_acc_rates = accretion_rates[np.abs(z_scores_accretion) < threshold]
    filtered_th_rates = theoretical_rates[np.abs(z_scores_theoretical) < threshold]

    # Also filter corresponding steps
    mask_acc = np.where(np.abs(z_scores_accretion) < threshold)[0]
    mask_th = np.where(np.abs(z_scores_theoretical) < threshold)[0]
    #print("\n", mask, "\n", len(filtered_acc_rates), " ", len(mask), len(filtered_th_rates))

    return filtered_acc_rates, filtered_th_rates, mask_acc, mask_th

# plots the mass of the star & mass accretion rates of the simulation
# input lists of masses, sound speed, timesteps, accretion rates, theoretical accretion rates & string to specify datatype
def plot_data(masses, sound_speeds, timesteps, accretion_rates, theoretical_rates, datatype, file_loc):
    plt.figure(figsize=(10, 10))

    # various x-axes for the plots & interpolation
    cumulative_time = np.cumsum(timesteps)
    spaced_timesteps = np.arange(0, 1000, 10)
    print("\n", len(cumulative_time))
    steps = np.arange(1, len(timesteps), 1)
    fine_time = np.linspace(min(cumulative_time), max(cumulative_time), 500)
    fine_steps = np.linspace(min(steps), max(steps), 100)

    # cubic spline interpolation 
    #cs_mass = CubicSpline(cumulative_time, masses)
    #cs_accretion = CubicSpline(cumulative_time[1:], accretion_rates)
    #cs_acc_steps = CubicSpline(steps, accretion_rates)
    #cs_th_steps = CubicSpline(steps, theoretical_rates)

    # univariate spline interpolation
    #f_acc_rates, f_th_rates, f_steps_acc, f_steps_th = filter_data(accretion_rates, theoretical_rates)
    #us_acc_steps = us(f_steps_acc, f_acc_rates, s=1)
    #us_th_steps = us(f_steps_th, f_th_rates, s=1)
    #us_mass = us(cumulative_time, masses, s=0.01)
    #us_accretion = us(cumulative_time[1:], accretion_rates, s=0.5)
    
    # plot mass evolution 
    plt.subplot(2, 2, 1)
    plt.scatter(spaced_timesteps, masses, label='Mass of Star', color='blue', marker='.')  # Cumulative sum of minDt
    #plt.plot(fine_time, cs_mass(fine_time), label='Cubic Spline Interpolation', color='violet', alpha=0.8)
    #plt.plot(fine_time, us_mass(fine_time), label='Univariate Interpolation', color='green', alpha=0.8)
    plt.title('Star Mass Over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Mass')
    plt.grid()
    plt.legend()

    # plot mass accretion rate, x-axis in realtime
    plt.subplot(2, 2, 2)
    plt.scatter(spaced_timesteps[1:], accretion_rates, label='Mass Accretion Rate', color='blue', marker='.')
    #plt.plot(fine_time[1:], cs_accretion(fine_time[1:]), label='Cubic Spline Interpolation', color='violet', alpha=0.8)
    #plt.plot(fine_time[1:], us_accretion(fine_time[1:]), label='Univariate Interpolation', color='green', alpha=0.8)
    plt.title('Mass Accretion Rate Over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Mass Accretion Rate')
    plt.ylim(0, 0.02)
    plt.grid()
    plt.legend()


    # plot a zoomed-in view of mass accretion rate, x-axis in realtime
    # plt.subplot(2, 2, 3)
    # plt.scatter(np.cumsum(timesteps)[1:], accretion_rates, label='Mass Accretion Rate', color='red', marker='.')
    # plt.title('Zoomed-In View of Mass Accretion Rate, (%s data)' %datatype)
    # plt.xlabel('Time (s)')
    # plt.ylabel('Mass Accretion Rate')
    # plt.ylim(min(accretion_rates), max(accretion_rates) * 0.04)
    # plt.grid()
    # plt.legend()

    # plot mass accretion rate compared to theoretical rate, x-axis in #timesteps
    plt.subplot(2, 2, 4)
    plt.scatter(np.arange(1, len(spaced_timesteps), 1), accretion_rates, label='Mass Accretion Rate (timesteps)', color='deeppink',marker='.' )
    # plt.plot(fine_steps, us_acc_steps(fine_steps), label='Univariate Interpolation (accretion rate)', color='green')
    plt.scatter(np.arange(1, len(spaced_timesteps), 1), theoretical_rates, label=f'Theoretical Mass Accretion Rate, alpha = {alpha}', color='blue', marker='.', alpha=0.8)
    #plt.plot(fine_steps, us_th_steps(fine_steps), label='Univariate Interpolation (theoretical rate)', color='green')
    plt.title('Mass Accretion Rate Compared to Analytical')
    plt.xlabel('timestep #')
    plt.ylabel('Mass Accretion Rate')
    plt.ylim(0, 0.02)
    plt.grid()
    plt.legend()

    fname = f'g_mass_accretion_{datatype}.png'
    plt.savefig(file_loc + fname)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    f5_m, f5_c = read_hdf5_data_2(f5, output, 10)
    txt_m, txt_c, timesteps = read_txt_data(file, 10)

    f5_acc = calc_mass_accretion(f5_m, timesteps)
    txt_acc = calc_mass_accretion(txt_m, timesteps)

    # print(len(f5_m), " ", len(timesteps), len(f5_acc))

    f5_th_acc = calc_theoretical_accretion(f5_c)
    #txt_th_acc = calc_theoretical_accretion(txt_c)

    plot_data(f5_m, f5_c, timesteps, f5_acc, f5_th_acc, "hdf5", plots)
    #plot_data(txt_m, txt_c, timesteps, txt_acc, txt_th_acc, "txt", plots)
