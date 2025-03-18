# script to plot mass accretion rate from simulation vs analytical mass accretion rate based on alpha * (c_s^3 / G)

import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import re
from scipy.optimize import minimize 

run_20 = './output/runs/run_disk_comb_20.hdf5'
run_20_beta = './output/runs/run_disk_comb_20_beta.hdf5'
run_20_mom = './output/runs/run_disk_mom_20.hdf5'
run_20_mom_beta = './output/runs/run_disk_mom_20_beta.hdf5'
run_radial = './output/runs/run_disk_radial_20.hdf5'
run_50_beta = './output/runs/run_disk_mom_50_beta.hdf5'
plots = './output/plots/mass-accretion/'
run_100_beta = '/home/lwatan/scratch/run_disk_mom_100_beta.hdf5'
run_100_comb_beta = '/home/lwatan/scratch/run_disk_comb_100_beta.hdf5'
run_100_radial_beta = '/home/lwatan/scratch/run_disk_radial_100_beta.hdf5'
run_mom_1e6 = '/home/lwatan/scratch/run_disk_mom_1e6_beta.hdf5'
run_comb_1e6 = '/home/lwatan/scratch/run_disk_comb_100_beta.hdf5'
run_rad_1e6 = '/home/lwatan/scratch/run_disk_radial_1e6_beta.hdf5'
run_comb_r1 = '/home/lwatan/scratch/run_disk_comb_1e6_r1_r1.hdf5'
run_half = '/home/lwatan/scratch/run_disk_half_1e6.hdf5'
run_double = '/home/lwatan/scratch/run_disk_double_1e6_short.hdf5'

# constants
G = 1.0 

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
        f5_time_list = time_data[::block_size]
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
    print("length of time array in calc_mass_acc: ", len(time_list))
    for i in range(1, len(time_list)):
        mass_rate = (masses[i] - masses[i - 1])/(time_list[i] - time_list[i-1])
        # print(time_list[i], " ", time_list[i]-time_list[i-1], "\n")
        # print(mass_rate, ",", time_steps[i], ",", mass_rate/time_steps[i], "\n", file=o)
        mass_accretion_rate.append(mass_rate)
    return mass_accretion_rate

def error_function(alpha, accretion_rate, sound_speeds): 
    analytical_rate =  alpha * (sound_speeds**3) / G 
    mse = np.mean((accretion_rate - analytical_rate)**2)
    return mse

def alpha_estimation(accretion_rates, sound_speeds, initial_alpha=0.4): 
    mse = lambda accretion_rates, analytical_rates: np.mean((accretion_rates - analytical_rates)**2)
    res = minimize(error_function, initial_alpha, args=(accretion_rates, sound_speeds), method='Nelder-Mead')
    best_alpha = res.x[0]
    return best_alpha
    
# plots the mass of the star & mass accretion rates of the simulation
# input lists of masses, sound speed, times, accretion rates, analytical accretion rates & string to specify datatype
def plot_data(masses, c, c_half, c_double, c_combined, time, accretion_rates, double_acc_rates, half_acc_rates, comb_acc_rates, datatype, file_loc):
    plt.figure(figsize=(10, 10))
    
    # plot mass evolution 
    plt.subplot(2, 2, 1)
    plt.scatter(time, masses, label='Mass of Star', color='blue', marker='.')  # Cumulative sum of minDt
    plt.title('Star Mass Over Time')
    plt.xlabel('Time (yr/2pi)')
    plt.ylabel('Mass (solar masses)')
    plt.grid()
    plt.legend()

    alpha = alpha_estimation(np.array(accretion_rates), np.array(c))
    analytical_rates = alpha * (c**3) / G
    # plot mass accretion rate compared to analytical rate, x-axis in #timesteps
    plt.subplot(2, 2, 2)
    plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate', color='blue',marker='.' )
    plt.scatter(time, analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha:.3f}', color='deeppink', marker='.', alpha=0.8)
    plt.title('Mass Accretion Rate Compared to Analytical')
    plt.xlabel(r'Time $\frac{{\left[ yr \right]}}{{\left[ 2\pi \right]}}$')
    plt.ylabel(r'Mass Accretion Rate $\left[ M_{{\odot}} \right]$')
    plt.ylim(0, 0.002)
    plt.grid()
    plt.legend()

    #alpha_comb = alpha_estimation(np.array(comb_acc_rates), np.array(c_combined))
    #an_rate_comb = alpha_comb * (c**3) / G
    # mass accretion rate of combined criteria compared to 1J accretion rate
    #plt.subplot(2, 2, 3)
    #plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate, momentum criterion', color='blue',marker='.' )
    #plt.scatter(time, analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha:.3f}', color='deeppink', marker='.', alpha=0.8)
    #plt.scatter(time[1:], comb_acc_rates, label=f'Mass Accretion Rate, distance & momentum criteria', color='dodgerblue', marker='.', alpha=0.8)
    #plt.scatter(time, an_rate_comb, label=f'Analytical Mass Accretion Rate of Combined Criteria, alpha = {alpha_comb:.3f}', color='orchid', marker='.', alpha=0.8)
    #plt.title('Mass Accretion Rates of single vs combined criteria')
    #plt.xlabel('Time (yr/2pi)')
    #plt.ylabel('Mass Accretion Rate (solar masses)')
    #plt.ylim(0, 0.002)
    #plt.grid()
    #plt.legend()

    # plot mass accretion rate compared to analytical rate, x-axis in #timesteps
    plt.subplot(2, 2, 4)
    plt.scatter(time[1:], double_acc_rates, label='Mass Accretion Rate (double threshold)', color='darkorange',marker='.' )
    plt.scatter(time[1:], half_acc_rates, label='Mass Accretion Rate (half threshold)', color='darkorchid',marker='.' )
    plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate', color='blue',marker='.', alpha=0.8 )
    plt.scatter(time[1:], analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha:.3f}', color='deeppink', marker='.', alpha=0.5)
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

def plot_mass_accretion(masses, accretion_rates, time, file_loc):
    plt.figure(figsize=(10, 10))
    
    # plot mass evolution 
    #plt.subplot(2, 1, 1)
    plt.scatter(time, masses, label='Mass of Star', color='blue', marker='.')  # Cumulative sum of minDt
    plt.title('Mass Over Time')
    plt.xlabel('Time (yr/2pi)')
    plt.ylabel('Mass (solar masses)')
    plt.grid()
    plt.legend()

    #plt.subplot(2, 1, 2)
    #plt.scatter(time[1:], accretion_rates, label='Mass Accretion Rate', color='blue',marker='.' )
    #plt.scatter(time[1:], analytical_rates, label=f'Analytical Mass Accretion Rate, alpha = {alpha}', color='deeppink', marker='.', alpha=0.8)
    #plt.title('Mass Accretion Rate Compared to Analytical')
    #plt.xlabel('Time (yr/2pi)')
    #plt.ylabel('Mass Accretion Rate (solar masses)')
    # plt.ylim(0, 0.002)
   # plt.grid()
    #plt.legend()

    fname = f'radial_accretion.png'
    plt.savefig(file_loc + fname)
    plt.tight_layout()
    plt.show()

def plot_fit(observed_rates, sound_speeds, best_alpha, time, file_loc, beta, run, G=1):
    analytical_rates = best_alpha * (sound_speeds**3) / G
    plt.figure(figsize=(10, 5))
    plt.scatter(time[1:], observed_rates, label="Modeled Accretion Rates", color="blue", alpha=0.7, marker='.')
    plt.plot(time, analytical_rates, label=f"Analytical Rates (alpha={best_alpha:.3f})", color="red", lw=2)
    plt.xlabel(r"Time $\frac{{\left[ yr \right]}}{{\left[ 2\pi \right]}}$")
    plt.ylabel(r"Mass Accretion Rate $\frac{{\left[ 2\pi \cdot M_{{\odot}} \right]}}{{\left[ yr \right]}}$")
    plt.ylim(0, 0.005)
    plt.legend()
    plt.grid()
    plt.title(f"Modeled vs Analytical Accretion Rates (beta={beta})")
    fname = f'analytical_accretion_{beta}_1e6_{run}_r1.pdf'
    plt.savefig(file_loc + fname)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    #comb_m_20, comb_c_20, comb_t_20 = read_hdf5_data_2(run_20_mom)
    #comb_m_20_beta, comb_c_20_beta, comb_t_20_beta = read_hdf5_data_2(run_20_mom_beta)
    #rad_m_20, rad_c_20, rad_t_20= read_hdf5_data_2(run_radial)
    #print("length of comb_20_t: ", len(comb_t_20))
    #print("length of comb_20_m: ", len(comb_m_20))
    #comb_acc_20 = calc_mass_accretion(comb_m_20, comb_t_20)
    #comb_acc_20_beta = calc_mass_accretion(comb_m_20_beta, comb_t_20_beta)
    #acc_rad = calc_mass_accretion(rad_m_20, rad_t_20)
    #best_alpha_20 = alpha_estimation(np.array(comb_acc_20), np.array(comb_c_20))
    #best_alpha_20_beta = alpha_estimation(np.array(comb_acc_20_beta), np.array(comb_c_20_beta))
    #best_alpha_rad = alpha_estimation(np.array(acc_rad),np.array(rad_c_20))

    #mom_m_50_beta, mom_c_50_beta, mom_t_50_beta = read_hdf5_data_2(run_50_beta)
    #mom_acc_50_beta = calc_mass_accretion(mom_m_50_beta, mom_t_50_beta)
    #best_alpha_50_beta = alpha_estimation(np.array(mom_acc_50_beta), np.array(mom_c_50_beta))

    runs = {
    #"run_comb_r1": run_comb_r1,
    #"run_half": run_half,
    "run_double": run_double
    }   
    for run_name, run_value in runs.items(): 
        run_type = run_name.split('_')[1]
        m, c, t = read_hdf5_data_2(run_value)
        acc = calc_mass_accretion(m, t)
        best_alpha = alpha_estimation(np.array(acc), np.array(c))

    # Plot the results
    #plot_fit(comb_acc_20, comb_c_20, best_alpha_20, comb_t_20, plots, "inf")
    #plot_fit(comb_acc_20_beta, comb_c_20_beta, best_alpha_20_beta, comb_t_20_beta, plots, "2pi")
        plot_fit(acc, c, best_alpha, t, plots, "2π", run_type)
    # plot_mass_accretion(rad_m_20, acc_rad, rad_t_20, plots)