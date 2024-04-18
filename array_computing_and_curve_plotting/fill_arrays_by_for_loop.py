import numpy as np
from array_computing_and_curve_plotting.fill_gaussian_lists import gaussian

# Initializing empty lists
x_values = []
y_values = []

if __name__ == '__main__':
    # Generating 41 uniformly spaced coordinates between -4 and 4
    x_values = np.linspace(-4, 4, 41)
    # Calculating Gaussian values for each x value
    y_values = [gaussian(x) for x in x_values]

    # Optionally, print out the values to verify
    for x, y in zip(x_values, y_values):
        print(f"x: {x:.1f}, y: {y:.6f}")
