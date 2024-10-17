import numpy as np
import matplotlib.pyplot as plt

file = '/Users/lilywatanabe/Desktop/eth/thesis/SPH-EXA-fork/output/star_mass.txt'

# Initialize lists to store mass and mass accretion rate
mass_list = []
time_steps = []

with open(file, 'r') as f:
    for line in f:
        # split the line by comma and strip whitespace
        time_step, mass_value = map(float, line.strip().split(','))
        time_steps.append(time_step)  # time difference per step (minDt)
        mass_list.append(mass_value)

# calculate mass accretion rate 
mass_accretion_rate = []
for i in range(1, len(mass_list)):
    mass_rate = (mass_list[i] - mass_list[i - 1]) / time_steps[i]  
    mass_accretion_rate.append(mass_rate)

# plottimng
time_steps_accretion = np.cumsum(time_steps)[1:]  # cumulative sum for real time

plt.figure(figsize=(12, 10))

# plot mass
plt.subplot(3, 1, 1)
plt.plot(np.cumsum(time_steps), mass_list, label='Mass of Star', color='blue')  # Cumulative sum of minDt
plt.title('Star Mass Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Mass (units)')
plt.grid()
plt.legend()

# plot mass accretion rate 
plt.subplot(3, 1, 2)
plt.plot(time_steps_accretion, mass_accretion_rate, label='Mass Accretion Rate', color='red')
plt.title('Mass Accretion Rate Over Time')
plt.xlabel('Time (s)')
plt.ylabel('Mass Accretion Rate')
plt.grid()
plt.legend()

# plot a zoomed-in view of mass accretion rate
plt.subplot(3, 1, 3)
plt.plot(time_steps_accretion, mass_accretion_rate, label='Mass Accretion Rate', color='red')
plt.title('Zoomed-In View of Mass Accretion Rate')
plt.xlabel('Time (s)')
plt.ylabel('Mass Accretion Rate')
plt.xlim(min(time_steps_accretion), max(time_steps_accretion))  
plt.ylim(min(mass_accretion_rate), max(mass_accretion_rate) * 0.1)  
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
