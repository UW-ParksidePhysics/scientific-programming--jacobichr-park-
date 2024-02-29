import numpy as np
import matplotlib.pyplot as plt

import math

pi = 3.1415926
angle = pi/4.
unit_value = math.sin(angle)**2 + math.cos(angle)**2
print(f'sin^2({angle:.3f}) + cos^2({angle:.3f}) = {unit_value:.3f}')