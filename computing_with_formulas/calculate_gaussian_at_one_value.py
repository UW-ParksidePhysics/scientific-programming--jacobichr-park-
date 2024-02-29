import numpy as np
import matplotlib.pyplot as plt
import math

m = 0 #mean
s = 2 #standard deviation
x = 1 #input_value
pi = math.pi
e = math.e

exponent = e**((((x - m)**2 / s**2)) * -0.5)
base = 1 / (s * math.sqrt(2*pi))
function = base * exponent

print(function)