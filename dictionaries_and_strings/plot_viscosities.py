import re
import numpy as np
import matplotlib.pyplot as plt

def parse_viscosity_data_adjusted(contents):
    """
    Parses the contents of the viscosity data file with an adjusted regex that accounts for
    multiple word gas names and returns a dictionary with the data.
    """
    viscosity_data_adjusted = {}
    # Adjusted regular expression to match lines with data including gas names with multiple words (This was something the AI reccomended while asking questions during bug testing)
    data_regex_adjusted = re.compile(
        r'^\s*(?P<gas>[a-zA-Z\s]+)\s+(?P<C>\d+)\s+(?P<T0>\d+\.\d+)\s+(?P<mu0>\d+\.\d+)\s*$'
    )

    for line in contents:
        match = data_regex_adjusted.match(line)
        if match:
            # Extract data for each gas
            gas_data = match.groupdict()
            gas_name = gas_data['gas'].strip()  # Remove any extra whitespace from the gas name
            # Convert strings to appropriate types
            gas_data['C'] = float(gas_data['C'])
            gas_data['T0'] = float(gas_data['T0'])
            gas_data['mu0'] = float(gas_data['mu0'])
            # Organize data in nested dictionary
            viscosity_data_adjusted[gas_name] = {
                'viscosity': {
                    'C': gas_data['C'],
                    'reference_temperature': gas_data['T0'],
                    'reference_viscosity': gas_data['mu0']
                }
            }
    return viscosity_data_adjusted

def calculate_viscosity(temperature, gas, viscosity_data):
    """
    Calculates the viscosity of a gas at a given temperature using the viscosity data.
    """
    # Retrieve the gas constants from the dictionary
    C = viscosity_data[gas]['viscosity']['C']
    T0 = viscosity_data[gas]['viscosity']['reference_temperature']
    mu0 = viscosity_data[gas]['viscosity']['reference_viscosity']

    # Apply the Sutherland's formula to calculate viscosity
    mu = mu0 * ((T0 + C) / (temperature + C)) * (temperature / T0)**1.5

    return mu

def plot_viscosities(temperature_range, gases, viscosity_data):
    """
    Plots the viscosity of given gases over a range of temperatures.
    """
    plt.figure(figsize=(10, 5))

    for gas in gases:
        # Calculate viscosity for each temperature in the range for the current gas
        viscosities = [calculate_viscosity(T, gas, viscosity_data) for T in temperature_range]
        plt.plot(temperature_range, viscosities, label=gas)

    plt.title('Viscosity of Gases over Temperature')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Viscosity ($10^{-6}$ Pa·s)')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
  # Read the content of the "viscosity_of_gases.dat" file
  with open('dictionaries_and_strings/viscosity_of_gases.dat', 'r') as file:
      contents = file.readlines()

  # Re-parse the data from the file contents using the adjusted function
  viscosity_data_adjusted = parse_viscosity_data_adjusted(contents)

  # Define the temperature range from 223 K to 373 K
  temperature_range = np.linspace(223, 373, 150)

  # Plot viscosities for air, carbon dioxide, and hydrogen
  plot_viscosities(temperature_range, ['air', 'carbon dioxide', 'hydrogen'], viscosity_data_adjusted)
