"""Generate Hertzsprung-Russell diagram like:
    https://en.wikipedia.org/wiki/Hertzsprung%E2%80%93Russell_diagram#/media/File:HRDiagram.png"""
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

astronomical_unit = 1.495978707e11  # meters
meters_to_light_years = 1./9.4607e15


def star_colormap(star_b_minus_vs):
    # Create color map from B-V = -0.33 (#7070ff) to 1.40 (#ff7f7f)
    # yellow = #ffff7f at B-V = 0.81
    number_of_gradient_points = 256
    white_index = int((0.33 / (0.33 + 1.40)) * number_of_gradient_points)
    yellow_index = int(((0.81 + .33) / (0.33 + 1.40)) * number_of_gradient_points)
    print(white_index, yellow_index, number_of_gradient_points)
    color_values = np.ones((number_of_gradient_points, 4))
    # Red values
    color_values[:white_index, 0] = np.linspace(112 / number_of_gradient_points,
                                                255 / number_of_gradient_points,
                                                white_index)
    color_values[white_index:yellow_index, 0] = np.linspace(255 / number_of_gradient_points,
                                                            255 / number_of_gradient_points,
                                                            yellow_index - white_index)
    color_values[yellow_index:, 0] = np.linspace(255 / number_of_gradient_points,
                                                 255 / number_of_gradient_points,
                                                 number_of_gradient_points - yellow_index)
    # Green values
    color_values[:white_index, 1] = np.linspace(112 / number_of_gradient_points,
                                                255 / number_of_gradient_points,
                                                white_index)
    color_values[white_index:yellow_index, 1] = np.linspace(255 / number_of_gradient_points,
                                                            255 / number_of_gradient_points,
                                                            yellow_index - white_index)
    color_values[yellow_index:, 1] = np.linspace(255 / number_of_gradient_points,
                                                 127 / number_of_gradient_points,
                                                 number_of_gradient_points - yellow_index)
    # Blue values
    color_values[:white_index, 2] = np.linspace(255 / number_of_gradient_points,
                                                255 / number_of_gradient_points,
                                                white_index)
    color_values[white_index:yellow_index, 2] = np.linspace(255 / number_of_gradient_points,
                                                            127 / number_of_gradient_points,
                                                            yellow_index - white_index)
    color_values[yellow_index:, 2] = np.linspace(127 / number_of_gradient_points,
                                                 127 / number_of_gradient_points,
                                                 number_of_gradient_points - yellow_index)
    new_colormap = ListedColormap(color_values)

    # Scale B-V values from 0 to 1
    scaled_b_minus_vs = (star_b_minus_vs - np.amin(star_b_minus_vs)) / (
            np.amax(star_b_minus_vs) - np.amin(star_b_minus_vs))

    return scaled_b_minus_vs, new_colormap


def parallax_to_distance(parallax):
    """Take parallax in milliarcseconds and convert to distance in meters"""
    parallax_in_radians = (parallax / 1000. / 3600.) * (2 * np.pi / 360.)
    distance = astronomical_unit / np.tan(parallax_in_radians)
    return distance


def apparent_to_absolute_magnitude(apparent_magnitude, distance):
    """Calculate absolute magnitude from apparent magnitude and distance in meters"""
    distance_in_parsecs = distance / (648000. * astronomical_unit / np.pi)
    absolute_magnitude = apparent_magnitude - 5*np.log10(distance_in_parsecs) + 5
    return absolute_magnitude


def read_hipparcos_data(file_lines):
  # Initialize an empty dictionary to store the data
  hipparcos_data = {}

  # Process each line
  for line in file_lines:
      # Remove newline characters and extra spaces, then split by space
      values = line.strip().split()

      # Process the values in groups of 4 (assuming each line is correctly formatted)
      for i in range(0, len(values), 4):
          # Extract the star data
          star_catalog_number = int(values[i])
          parallax = float(values[i + 1])
          apparent_magnitude = float(values[i + 2])
          blue_minus_visual = float(values[i + 3])

          # Add the star data to the dictionary
          hipparcos_data[star_catalog_number] = {
            'parallax': parallax,
            'apparent_magnitude': apparent_magnitude,
            'blue_minus_visual': blue_minus_visual
        }

  return hipparcos_data

if __name__ == '__main__':
  # Read the content of the "hipparcos_data.txt" file
  with open('dictionaries_and_strings/hipparcos_data.txt', 'r') as file:
      hipparcos_data_lines = file.readlines()

  # Parse the hipparcos data using the defined function
  hipparcos_data_parsed = read_hipparcos_data(hipparcos_data_lines)

  # Convert parsed data into NumPy arrays for processing
  parallaxes = np.array([data['parallax'] for data in hipparcos_data_parsed.values()])
  apparent_magnitudes = np.array([data['apparent_magnitude'] for data in hipparcos_data_parsed.values()])
  blue_minus_visuals = np.array([data['blue_minus_visual'] for data in hipparcos_data_parsed.values()])

  # Calculate distances and absolute magnitudes
  distances = parallax_to_distance(parallaxes)
  absolute_magnitudes = apparent_to_absolute_magnitude(apparent_magnitudes, distances)

  # Use dark style for plot
  plt.style.use('dark_background')

  # Reverse the absolute magnitude so that negative values appear on top
  reversed_absolute_magnitudes = np.negative(absolute_magnitudes)

  # Get color map to match star colors
  scaled_b_minus_v, hr_diagram_colormap = star_colormap(blue_minus_visuals)

  # Create axes labels and make the axes identical to the model figure
  plt.xlabel('B-V Color Index')
  plt.ylabel('Absolute Magnitude')
  plt.xlim(-0.5, 2.5)  # Assuming a range from the model image
  plt.ylim(15, -10)    # Assuming a range from the model image

  # Scatter plot of B-V vs absolute magnitude with the defined scatter marker size
  plt.scatter(blue_minus_visuals, reversed_absolute_magnitudes, c=scaled_b_minus_v, cmap=hr_diagram_colormap, s=1)

  # Add creator attribution text
  plt.text(0.05, 0.05, 'Created by Jacob, Noah', transform=plt.gca().transAxes, fontsize=8, color='white')

  plt.show()