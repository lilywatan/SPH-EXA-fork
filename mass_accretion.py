# script to plot mass accretion rate from simulation vs analytical mass accretion rate based on alpha * (c_s^3 / G)

import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import re

f5 = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_1J_5000.hdf5'
f5_double = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_2J_5000.hdf5'
f5_half = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_05J_5000.hdf5'
combined_criteria = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_disk_1J_5000_dist.hdf5'
output = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/5000.txt'
plots = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/plots/'
cloud = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/run_cloud_1000.hdf5'

# constants
G = 1 
alpha = 0.4 # scaling parameter for analytical mass accretion rate

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

        # reshape arrays into chunks for averaging (mean of n datapoints)
        f5_mass_list = np.mean(star_mass_data[:num_chunks * block_size].reshape(-1, block_size), axis=1)
        f5_sound_speed = np.mean(sound_speed_data[:num_chunks * block_size].reshape(-1, block_size), axis=1)
        # list of every n-th time 
        f5_time_list = time_data[:-1:block_size]
        # print("length of mass_list: ", len(f5_mass_list), " ")
        # print("length of time_list: ", len(f5_time_list), " ")

        # write to output if specified
        if output is not None:
            with open(output, 'a') as o:
                for i in range(1, len(time_data), 1):
                    o.write(f"{time_data[i]}, {time_data[i] - time_data[i-1]} \n")

                o.write("averaged times (every 10th timestep): \n")
                for t in zip(f5_time_list):
                    o.write(f"{t} \n")

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
        analytical_rate.append(alpha * ((c_s[i]**3) / G)) # M_dot = a * c_s^3/G
    return analytical_rate

# plots the mass of the star & mass accretion rates of the simulation
# input lists of masses, sound speed, times, accretion rates, analytical accretion rates & string to specify datatype
def plot_data(masses, sound_speeds, time, accretion_rates, analytical_rates, double_acc_rates, half_acc_rates, comb_acc_rates, datatype, file_loc):
    plt.figure(figsize=(10, 10))
    
    # plot mass evolution 
    plt.subplot(2, 2, 1)
    plt.scatter(time, masses, label='Mass of Star', color='blue', marker='.')  # Cumulative sum of minDt
    plt.title('Star Mass Over Time')
    plt.xlabel('Time (yr/2pi)')
    plt.ylabel('Mass (solar masses)')
    plt.grid()
    plt.legend()

    # plot mass accretion rate, x-axis in realtime
    #plt.subplot(2, 2, 2)
    #plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate', color='blue', marker='.')
    #plt.title('Mass Accretion Rate Over Time')
    #plt.xlabel('Time (yr)')
    #plt.ylabel('Mass Accretion Rate (solar masses)')
    #plt.ylim(0, 0.002)
    #plt.grid()
    #plt.legend()

    # plot mass accretion rate compared to analytical rate, x-axis in #timesteps
    plt.subplot(2, 2, 2)
    plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate', color='blue',marker='.' )
    plt.scatter(time[1:], analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha}', color='deeppink', marker='.', alpha=0.8)
    plt.title('Mass Accretion Rate Compared to Analytical')
    plt.xlabel('Time (yr/2pi)')
    plt.ylabel('Mass Accretion Rate (solar masses)')
    plt.ylim(0, 0.002)
    plt.grid()
    plt.legend()

    # mass accretion rate of combined criteria compared to 1J accretion rate
    plt.subplot(2, 2, 3)
    plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate, momentum criterion', color='blue',marker='.' )
    plt.scatter(time[1:], comb_acc_rates, label=f'Mass Accretion Rate, distance & momentum criteria', color='deeppink', marker='.', alpha=0.8)
    plt.title('Mass Accretion Rates of single vs combined criteria')
    plt.xlabel('Time (yr/2pi)')
    plt.ylabel('Mass Accretion Rate (solar masses)')
    plt.ylim(0, 0.002)
    plt.grid()
    plt.legend()

    # plot mass accretion rate compared to analytical rate, x-axis in #timesteps
    plt.subplot(2, 2, 4)
    plt.scatter(time[1:], double_acc_rates, label='Mass Accretion Rate (double threshold)', color='darkorange',marker='.' )
    plt.scatter(time[1:], half_acc_rates, label='Mass Accretion Rate (half threshold)', color='darkorchid',marker='.' )
    plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate', color='blue',marker='.', alpha=0.8 )
    plt.scatter(time[1:], analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha}', color='deeppink', marker='.', alpha=0.5)
    plt.title('Mass Accretion Rate Compared to Analytical')
    plt.xlabel('Time (yr/2pi)')
    plt.ylabel('Mass Accretion Rate (solar masses)')
    plt.ylim(0, 0.002)
    plt.grid()
    plt.legend()

    fname = f'mass_accretion_multJ.png'
    plt.savefig(file_loc + fname)
    plt.tight_layout()
    plt.show()

def plot_mass_accretion(masses, accretion_rates, analytical_rates, time, file_loc):
    plt.figure(figsize=(10, 10))
    
    # plot mass evolution 
    plt.subplot(2, 1, 1)
    plt.scatter(time, masses, label='Mass of Star', color='blue', marker='.')  # Cumulative sum of minDt
    plt.title('Mass Over Time')
    plt.xlabel('Time (yr/2pi)')
    plt.ylabel('Mass (solar masses)')
    plt.grid()
    plt.legend()

    plt.subplot(2, 1, 2)
    plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate', color='blue',marker='.' )
    plt.scatter(time[1:], analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha}', color='deeppink', marker='.', alpha=0.8)
    plt.title('Mass Accretion Rate Compared to Analytical')
    plt.xlabel('Time (yr/2pi)')
    plt.ylabel('Mass Accretion Rate (solar masses)')
    # plt.ylim(0, 0.002)
    plt.grid()
    plt.legend()

    fname = f'cloud_accretion.png'
    plt.savefig(file_loc + fname)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    f5_m, f5_c, f5_t = read_hdf5_data_2(f5, output, 10)
    double_m, double_c, double_t = read_hdf5_data_2(f5_double, None, 10)
    half_m, half_c, half_t = read_hdf5_data_2(f5_half, None, 10)
    cloud_m, cloud_c, cloud_t = read_hdf5_data_2(cloud, output, 10)
    comb_m, comb_c, comb_t = read_hdf5_data_2(combined_criteria, output, 10)

    f5_acc = calc_mass_accretion(f5_m, f5_t)
    double_acc = calc_mass_accretion(double_m, double_t)
    half_acc = calc_mass_accretion(half_m, half_t)
    cloud_acc = calc_mass_accretion(cloud_m, cloud_t)
    combined_acc = calc_mass_accretion(comb_m, comb_t)

    f5_an_acc = calc_analytical_accretion(f5_c)
    double_an_acc = calc_analytical_accretion(double_c)
    half_an_acc = calc_analytical_accretion(half_c)
    cloud_an_acc = calc_analytical_accretion(cloud_c)

    plot_data(f5_m, f5_c, f5_t, f5_acc, f5_an_acc, double_acc, half_acc, combined_acc, "hdf5", plots)
    plot_mass_accretion(cloud_m, cloud_acc, cloud_an_acc, cloud_t, plots)