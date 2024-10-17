import numpy as np
import matplotlib.pyplot as plt
import sys

file = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/star_mass_sound_speed.txt'
output = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/data_analysis_theoretical.txt'

# constants
G = 6.6726e-8 # taken from one of the sphexa files (cooling?)
alpha = 0.1 # scaling parameter for theoretical mass accretion rate

# Initialize lists to store mass and mass accretion rate
mass_list = []
time_steps = []
sound_speed = []

with open(output, 'a') as o:
    with open(file, 'r') as f:
        for line in f:
            # split the line by comma and strip whitespace
            c_s, time_step, mass_value = map(float, line.strip().split(','))
            sound_speed.append(c_s)
            time_steps.append(time_step)  # time difference per step (minDt)
            mass_list.append(mass_value)

    # calculate mass accretion rate
    mass_accretion_rate = []
    for i in range(1, len(mass_list)):
        mass_rate = (mass_list[i] - mass_list[i - 1]) / time_steps[i]
        print(mass_rate, ",", time_steps[i], ",", mass_rate/time_steps[i], "\n", file=o)
        mass_accretion_rate.append(mass_rate)

    # calculate theoretical mass accretion rate
    theoretical_rate = []
    for i in range(1, len(sound_speed)):
        theoretical_rate.append(G / (alpha * sound_speed[i] * sound_speed[i] * sound_speed[i]))

# plot timing
time_steps_accretion = np.cumsum(time_steps)[1:]  # cumulative sum for real time

plt.figure(figsize=(12, 10))

# plot mass
plt.subplot(2, 2, 1)
plt.scatter(np.cumsum(time_steps), mass_list, label='Mass of Star', color='blue', marker='.')  # Cumulative sum of minDt
plt.title('Star Mass Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Mass (units)')
plt.grid()
plt.legend()

# plot mass accretion rate with timesteps
plt.subplot(2, 2, 2)
plt.scatter(time_steps_accretion, mass_accretion_rate, label='Mass Accretion Rate', color='red', marker='.')
plt.title('Mass Accretion Rate Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Mass Accretion Rate')
plt.grid()
plt.legend()

# plot a zoomed-in view of mass accretion rate
plt.subplot(2, 2, 3)
plt.scatter(time_steps_accretion, mass_accretion_rate, label='Mass Accretion Rate', color='red', marker='.')
plt.title('Zoomed-In View of Mass Accretion Rate')
plt.xlabel('Time (s)')
plt.ylabel('Mass Accretion Rate')
plt.xlim(min(time_steps_accretion), max(time_steps_accretion))  
plt.ylim(min(mass_accretion_rate), max(mass_accretion_rate) * 0.1)
plt.grid()
plt.legend()

# plot mass accretion rate without timesteps
plt.subplot(2, 2, 4)
plt.scatter(np.arange(1, len(time_steps), 1), mass_accretion_rate, label='Mass Accretion Rate no timesteps', color='violet',marker='.' )
plt.scatter(np.arange(1, len(time_steps), 1), theoretical_rate, label='Theoretical Mass Accretion Rate', color='blue', marker='.')
plt.title('Mass Accretion Rate Over Time')
plt.xlabel('timestep #')
plt.ylabel('Mass Accretion Rate')
plt.ylim(0, 0.02)
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
