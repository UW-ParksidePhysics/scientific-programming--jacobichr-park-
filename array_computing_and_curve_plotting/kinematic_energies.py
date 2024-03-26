# define kinetic energy function
# define potential energy function
# Import height annd velocity functions of time from freefall_kinematics.py
# define the array of times from 0 to 2 v_0/g
# plot the energies versus time
# predefine mass and inital velocity
import numpy as np
import matplotlib.pyplot as plt
from freefall_kinematics import height, velocity


def kinetic_energy(velocity, mass):
  return 0.5 * mass * velocity**2


def potential_energy(height, mass, gravitational_acceleration):
  return mass * gravity * height


def calculate_energies(mass, initial_velocity, gravitational_acceleration=9.8):
  end_time = 2 * initial_velocity / gravitational_acceleration
  times = np.linspace(0, end_time, 100)
  heights = height(times, initial_height, initial_velocity, gravitational_acceleration)
  velocities = velocity(times, initial_velocity, gravitational_acceleration)
  kinetic_energies = kinetic_energy(velocities, mass)
  potential_energies = potential_energy(heights, mass, gravitational_acceleration=gravitational_acceleration)
  total_energies = kinetic_energies + potential_energies
  return [times, kinetic_energies, potential_energies, total_energies]


def plot_energies(times, kinetic_energies, potential_energies, total_energies):
  plt.plot(times, kinetic_energies, label=r'$K(t)$')
  plt.plot(times, potential_energies, label=r'$P(t)$')
  plt.plot(times, total_energies, label=r'$E(t)$')
  plt.xlabel('$t$ (s)')
  plt.ylabel('$E, K, P$ (J)')
  plt.show()
  return


if __name__ == "__main__":
  object_mass = 1.
  starting_height = 1.
  starting_velocity = 1.
  time_sample_number = 100
  calculate_energies(object_mass, starting_height, starting_velocity, time_sample_number)