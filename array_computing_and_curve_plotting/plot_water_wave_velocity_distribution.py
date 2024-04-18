import numpy as np
import matplotlib.pyplot as plt

# Constants
g = 9.81  # Acceleration due to gravity (m/s^2)
s_pure = 7.2e-2  # Surface tension of pure water (N/m)
s_sea = 6.5e-2  # Surface tension of seawater (N/m)
p_pure = 998  # Density of pure water (kg/m^3)
p_sea = 1025  # Density of seawater (kg/m^3)
h = 50  # Water depth (m)

# Wave speed function
def c_lambda(l, s, p):
    return np.sqrt((g*l/(2*np.pi)) * (1 + s*(4*np.pi**2)/(p*g*l)) * np.tanh((2*np.pi*h)/l))


if __name__ == "__main__":
  # Wavelength ranges
  small_lambda_range = np.linspace(0.001, 0.1, 500)
  large_lambda_range = np.linspace(1, 2000, 500)

  # Calculate wave speeds
  c_small_pure = c_lambda(small_lambda_range, s_pure, p_pure)
  c_small_sea = c_lambda(small_lambda_range, s_sea, p_sea)
  c_large_pure = c_lambda(large_lambda_range, s_pure, p_pure)
  c_large_sea = c_lambda(large_lambda_range, s_sea, p_sea)

  # Plotting
  fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 6))

  # Small λ plot
  ax1.plot(small_lambda_range, c_small_pure, label='Pure Water')
  ax1.plot(small_lambda_range, c_small_sea, label='Seawater')
  ax1.set_title('Wave Speed c(λ) for Small λ')
  ax1.set_xlabel('Wavelength λ (m)')
  ax1.set_ylabel('Wave Speed c(λ) (m/s)')
  ax1.legend()
  ax1.grid(True)

  # Large λ plot
  ax2.plot(large_lambda_range, c_large_pure, label='Pure Water')
  ax2.plot(large_lambda_range, c_large_sea, label='Seawater')
  ax2.set_title('Wave Speed c(λ) for Large λ')
  ax2.set_xlabel('Wavelength λ (m)')
  ax2.set_ylabel('Wave Speed c(λ) (m/s)')
  ax2.set_xscale('log')  # Log scale for better visualization of large range
  ax2.legend()
  ax2.grid(True)

  plt.tight_layout()
  plt.show()
