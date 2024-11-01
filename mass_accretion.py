# script to plot mass accretion rate from simulation vs analytical mass accretion rate based on alpha * (c_s^3 / G)

import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
from scipy.interpolate import CubicSpline 
from scipy.interpolate import UnivariateSpline as us
from scipy.stats import zscore
import re

f5 = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_1J_5000.hdf5'
f5_double = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_2J_5000.hdf5'
f5_half = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_05J_5000.hdf5'
output = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/5000.txt'
plots = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/plots/'

# constants
G = 1 
G2 = G * G
alpha = 0.5 # scaling parameter for analytical mass accretion rate

def read_hdf5_data_2(file, output=None, block_size=1):
    f5_mass_list = []
    f5_sound_speed = []

    with h5py.File(file, 'r') as f: 
        step_keys = sorted(f.keys(), key=lambda x: int(re.search(r'\d+', x).group()))

        star_mass_data = []
        sound_speed_data = []
        time = []
        
        # extract all star mass and sound speed data
        for step_key in step_keys:
            star_mass_data.append(f[step_key].attrs['star::m'])
            sound_speed_data.append(np.mean(f[step_key]['c']))
            time.append(f[step_key].attrs['time'])

        # convert lists to numpy arrays for vectorized operations
        star_mass_data = np.array(star_mass_data)
        sound_speed_data = np.array(sound_speed_data)
        time_data = np.array(time)
        print("length of time_data: ", len(time_data), " ")

        # determine number of full chunks
        num_chunks = len(star_mass_data) // block_size

        # reshape arrays into chunks for averaging
        f5_mass_list = np.mean(star_mass_data[:num_chunks * block_size].reshape(-1, block_size), axis=1)
        f5_sound_speed = np.mean(sound_speed_data[:num_chunks * block_size].reshape(-1, block_size), axis=1)
        f5_time_list = np.mean(time_data[:num_chunks * block_size].reshape(-1, block_size), axis=1)

        # write to output if specified
        if output is not None:
            with open(output, 'a') as o:
                #for mass, speed in zip(f5_mass_list, f5_sound_speed):
                    #o.write(f"{mass}, {speed} \n")
                for i in range(1, len(time_data), 1):
                    o.write(f"{time_data[i]}, {time_data[i] - time_data[i-1]} \n")

                o.write("averaged times (every 10th timestep): \n")
                for t in zip(f5_time_list):
                    o.write(f"{t} \n")

        # print(len(f5_mass_list), " ", len(f5_sound_speed))

    return f5_mass_list, f5_sound_speed, f5_time_list


# calculate mass accretion rate from simulation data
# input: lists of masses & timesteps, return list of mass accretion rates
def calc_mass_accretion(masses, time_list):
    mass_accretion_rate = []
    for i in range(1, len(masses)):
        mass_rate = (masses[i] - masses[i - 1])/(time_list[i] - time_list[i-1])
        # print(time_list[i], " ", time_list[i]-time_list[i-1], "\n")
        # print(mass_rate, ",", time_steps[i], ",", mass_rate/time_steps[i], "\n", file=o)
        mass_accretion_rate.append(mass_rate)
    return mass_accretion_rate

# calculates the analytical mass accretion rate
# input: list of sound speeds, return: list of analytical accretion values
def calc_analytical_accretion(c_s):
    analytical_rate = []
    for i in range(1, len(c_s)):
        analytical_rate.append(alpha * ((c_s[i] * c_s[i] * c_s[i]) / G)) # M_dot = a * c_s^3/G
    return analytical_rate

def filter_data(accretion_rates, analytical_rates):
    accretion_rates = np.array(accretion_rates)
    analytical_rates = np.array(analytical_rates)

    # Compute Z-scores for accretion rates
    z_scores_accretion = zscore(accretion_rates)
    z_scores_analytical = zscore(analytical_rates)

    # Define a threshold for Z-score to detect outliers (e.g., 3 standard deviations)
    threshold = 2

    # Filter out outliers
    filtered_acc_rates = accretion_rates[np.abs(z_scores_accretion) < threshold]
    filtered_an_rates = analytical_rates[np.abs(z_scores_analytical) < threshold]

    # Also filter corresponding steps
    mask_acc = np.where(np.abs(z_scores_accretion) < threshold)[0]
    mask_an = np.where(np.abs(z_scores_analytical) < threshold)[0]
    #print("\n", mask, "\n", len(filtered_acc_rates), " ", len(mask), len(filtered_an_rates))

    return filtered_acc_rates, filtered_an_rates, mask_acc, mask_an

# plots the mass of the star & mass accretion rates of the simulation
# input lists of masses, sound speed, timesteps, accretion rates, analytical accretion rates & string to specify datatype
def plot_data(masses, sound_speeds, timesteps, accretion_rates, analytical_rates, double_acc_rates, half_acc_rates, datatype, file_loc):
    plt.figure(figsize=(10, 10))

    # various x-axes for the plots & interpolation
    cumulative_time = np.cumsum(timesteps)
    spaced_timesteps = np.arange(0, 1000, 10)
    # print("\n", len(cumulative_time))
    steps = np.arange(1, len(timesteps), 1)
    fine_time = np.linspace(min(cumulative_time), max(cumulative_time), 500)
    fine_steps = np.linspace(min(steps), max(steps), 100)
    
    # plot mass evolution 
    plt.subplot(2, 2, 1)
    plt.scatter(timesteps, masses, label='Mass of Star', color='blue', marker='.')  # Cumulative sum of minDt
    #plt.plot(fine_time, cs_mass(fine_time), label='Cubic Spline Interpolation', color='violet', alpha=0.8)
    #plt.plot(fine_time, us_mass(fine_time), label='Univariate Interpolation', color='green', alpha=0.8)
    plt.title('Star Mass Over Time')
    plt.xlabel('Time')
    plt.ylabel('Mass')
    plt.grid()
    plt.legend()

    # plot mass accretion rate, x-axis in realtime
    plt.subplot(2, 2, 2)
    plt.scatter(timesteps[1:], accretion_rates, label='Mass Accretion Rate', color='blue', marker='.')
    #plt.plot(fine_time[1:], cs_accretion(fine_time[1:]), label='Cubic Spline Interpolation', color='violet', alpha=0.8)
    #plt.plot(fine_time[1:], us_accretion(fine_time[1:]), label='Univariate Interpolation', color='green', alpha=0.8)
    plt.title('Mass Accretion Rate Over Time')
    plt.xlabel('Time')
    plt.ylabel('Mass Accretion Rate')
    plt.ylim(0, 0.002)
    plt.grid()
    plt.legend()

    # plot mass accretion rate compared to analytical rate, x-axis in #timesteps
    plt.subplot(2, 2, 3)
    plt.scatter(timesteps[1:], accretion_rates, label='Mass Accretion Rate', color='blue',marker='.' )
    # plt.plot(fine_steps, us_acc_steps(fine_steps), label='Univariate Interpolation (accretion rate)', color='green')
    plt.scatter(timesteps[1:], analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha}', color='deeppink', marker='.', alpha=0.8)
    #plt.plot(fine_steps, us_an_steps(fine_steps), label='Univariate Interpolation (analytical rate)', color='green')
    plt.title('Mass Accretion Rate Compared to Analytical')
    plt.xlabel('Time')
    plt.ylabel('Mass Accretion Rate')
    plt.ylim(0, 0.002)
    plt.grid()
    plt.legend()

    # plot mass accretion rate compared to analytical rate, x-axis in #timesteps
    plt.subplot(2, 2, 4)
    plt.scatter(timesteps[1:], double_acc_rates, label='Mass Accretion Rate (double threshold)', color='darkorange',marker='.' )
    plt.scatter(timesteps[1:], half_acc_rates, label='Mass Accretion Rate (half threshold)', color='darkorchid',marker='.' )
    plt.scatter(timesteps[1:], accretion_rates, label='Mass Accretion Rate', color='blue',marker='.', alpha=0.8 )
    # plt.plot(fine_steps, us_acc_steps(fine_steps), label='Univariate Interpolation (accretion rate)', color='green')
    plt.scatter(timesteps[1:], analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha}', color='deeppink', marker='.', alpha=0.5)
    #plt.plot(fine_steps, us_an_steps(fine_steps), label='Univariate Interpolation (analytical rate)', color='green')
    plt.title('Mass Accretion Rate Compared to Analytical')
    plt.xlabel('Time')
    plt.ylabel('Mass Accretion Rate')
    plt.ylim(0, 0.002)
    plt.grid()
    plt.legend()

    fname = f'g_mass_accretion_{datatype}_multJ.png'
    plt.savefig(file_loc + fname)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    f5_m, f5_c, f5_t = read_hdf5_data_2(f5, output, 10)
    double_m, double_c, double_t = read_hdf5_data_2(f5_double, None, 10)
    half_m, half_c, half_t = read_hdf5_data_2(f5_half, None, 10)

    f5_acc = calc_mass_accretion(f5_m, f5_t)
    double_acc = calc_mass_accretion(double_m, double_t)
    half_acc = calc_mass_accretion(half_m, half_t)

    f5_an_acc = calc_analytical_accretion(f5_c)
    double_an_acc = calc_analytical_accretion(double_c)
    half_an_acc = calc_analytical_accretion(half_c)

    plot_data(f5_m, f5_c, f5_t, f5_acc, f5_an_acc, double_acc, half_acc, "hdf5", plots)
