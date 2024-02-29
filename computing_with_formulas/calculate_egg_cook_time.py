import numpy as np
import matplotlib.pyplot as plt
import math

mass_small = 47 #grams, mass of small egg
mass_large = 67 #grams, mass of large egg
density = 1.038 #grams/cm^3, denisty of both eggs
specific_heat = 3.7 #J/g*C, specific heat of both eggs
thermal_conductivity = 5.4 * 10**(-3) #W/cm*K, thermal conductivity of both eggs

Tw = 100 #C, temperature of boiling water
Ty = 70 #C, temperature of boiling yolk
Tf = 4 #C, temperature of fridge
Tr = 20 #C, temperature of room

naturallog_room = math.log(0.76 * ((Tr-Tw) / (Ty-Tw)))
naturallog_fridge = math.log(0.76 * ((Tf-Tw) / (Ty-Tw)))
base_small = (mass_small**(2/3) * specific_heat * density**(1/3)) / (thermal_conductivity * math.pi**2 * (4*math.pi / 3)**(2/3))
base_large = (mass_large**(2/3) * specific_heat * density**(1/3)) / (thermal_conductivity * math.pi**2 * (4*math.pi / 3)**(2/3))

print(f'Time small egg = {base_small * naturallog_room:.2f} s')
print(f'Time large egg = {base_large * naturallog_room:.2f} s')