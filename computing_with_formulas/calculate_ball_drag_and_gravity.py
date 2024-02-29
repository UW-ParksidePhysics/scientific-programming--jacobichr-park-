import numpy as np
import matplotlib.pyplot as plt
import math

pi = math.pi

air_density = 1.2 #kg/m^3, air density
ball_radius = 0.11 #meters, radius of ball
ball_mass = 0.43 #kg, mass of ball
drag_coefficient = 0.2 #drag coefficient
ball_velocity_1 = 33.33333 #m/s, initial velocity hard kick
ball_velocity_2 = 2.777777 #m/s, initial velocity soft kick
gravitational_acceleration = 9.81 #m/s^2, gravity

cross_area = math.pi * (ball_radius**2) # Cross sectional area #
gravitational_force = ball_mass * gravitational_acceleration
drag_force_1 = 0.5 * air_density * cross_area * drag_coefficient * ball_velocity_1**2 # drag force on hard kick #
drag_force_2 = 0.5 * air_density * cross_area * drag_coefficient * ball_velocity_2**2 # drag force on soft kick #

print(f'cross area = {cross_area:.2f} m^2')
print(f'gravitational force = {gravitational_force:.2f} N')
print(f'Drag force hard kick = {drag_force_1}N')
print(f'Drag force soft kick = {drag_force_2}N')