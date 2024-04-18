import numpy as np
from array_computing_and_curve_plotting.fill_gaussian_lists import gaussian

# Initializing empty lists
x_values = []
y_values = []

if __name__ == '__main__':
  # Generating 41 uniformly spaced coordinates between -4 and 4 as a NumPy array
  # Note: I had already been using np.linespace as it was suggested by the AI.
  x_values = np.linspace(-4, 4, 41)
  # Evaluating Gaussian values for the entire array of x values
  y_values = gaussian(x_values)

  # Optionally, print out the values to verify
  for x, y in zip(x_values, y_values):
      print(f"x: {x:.1f}, g(x): {y:.6f}")
