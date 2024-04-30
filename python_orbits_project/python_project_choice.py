from skyfield.api import load, Loader
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

# Initialize the loader and timescale
loader = Loader('~/.skyfield-data', verbose=False)
ts = loader.timescale()

def get_date_input(prompt, ts):
    while True:
        print(prompt)
        try:
            year = int(input("Enter the year, must be from 1550-2650 or 1969-17191 (AD): "))
            month = int(input("Enter the month (MM), e.g., 07 for July: "))
            day = int(input("Enter the day (DD), e.g., 15: "))
            # For BC dates, we need to subtract one because there is no year 0 in the astronomical year numbering.
            year = year if year > 0 else year - 1
            date = ts.utc(year, month, day)
            return date
        except ValueError as e:
            print(f"Invalid date entered: {e}")

def select_ephemeris_file(date1, date2, ephemeris_files):
    # Use .utc[0] to extract the year for the comparison
    year1, year2 = date1.utc[0], date2.utc[0]
    for file_name, (start_year, end_year) in ephemeris_files.items():
        if start_year <= year1 <= end_year and start_year <= year2 <= end_year:
            return file_name
    return None

def calculate_celestial_distance(ephemeris, celestial_body_name, date):
    sun = ephemeris['sun']
    celestial_body = ephemeris[celestial_body_name]
    distance = sun.at(date).observe(celestial_body).distance().au
    return distance

def calculate_velocities(dates, distances):
  velocities = [None]  # No velocity for the first date
  for i in range(1, len(dates) - 1):
      # Central difference method: (f(x+h) - f(x-h)) / (2h)
      time_diff = (dates[i + 1].tt - dates[i - 1].tt) * 24 * 3600  # Time difference in seconds
      distance_diff = distances[i + 1] - distances[i - 1]  # Distance difference in AU
      velocity = distance_diff / time_diff  # Velocity in AU/s
      velocities.append(velocity * 149597870.7 * 1000)  # Convert AU/s to m/s
  velocities.append(None)  # No velocity for the last date
  return velocities

def calculate_gravitational_force(mass_body, mass_sun, distance):
  G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
  force = G * mass_body * mass_sun / (distance * 149597870.7e3)**2  # Convert AU to meters for the distance
  return force

def calculate_acceleration(force, mass_body):
  acceleration = (force / mass_body) * 1000000  # Acceleration in micrometers/s^2
  return acceleration


def print_distance_table(years, neptune_distances, neptune_velocities, neptune_accelerations, pluto_distances, pluto_velocities, pluto_accelerations):
  table = zip(years, neptune_distances, neptune_velocities, neptune_accelerations, pluto_distances, pluto_velocities, pluto_accelerations)
  headers = ["Year", "Neptune Distance (AU)", "Neptune RoC (m/s)", "Neptune Acceleration (µ/s²)", "Pluto Distance (AU)", "Pluto RoC (m/s)", "Pluto Acceleration (µ/s²)"]
  print(tabulate(table, headers=headers, floatfmt=".2f", missingval="N/A"))

def print_closer_intervals(year_labels, neptune_distances, pluto_distances):
  intervals = []
  i = 1
  while i < len(year_labels) - 1:
      # Find the start of an interval where Pluto is closer
      if pluto_distances[i] < neptune_distances[i]:
          start_year = year_labels[i]
          while i < len(year_labels) - 1 and pluto_distances[i] < neptune_distances[i]:
              i += 1
          end_year = year_labels[i]
          intervals.append((start_year, end_year))
      i += 1

  print("\nIntervals when Pluto was closer to the Sun than Neptune:")
  for start, end in intervals:
      print(f"From {start:.2f} to {end:.2f}")



# Define ephemeris files with their supported year ranges
ephemeris_files = {
    'de440s.bsp': (1849, 2150),
    'de440.bsp': (1550, 2650),
    'de441.bsp': (1969, 17191)
}

date1 = get_date_input("Enter the first date:", ts)
date2 = get_date_input("Enter the second date:", ts)

# Ensure date1 is before date2
if date2.tt < date1.tt:
    date1, date2 = date2, date1

file_name = select_ephemeris_file(date1, date2, ephemeris_files)
if file_name:
    ephemeris = loader(file_name)
    print(f"Loaded {file_name} for the date range.")
else:
    print("No ephemeris file covers the entire date range entered.")
    exit()

# Calculate intervals and dates for plotting
days_difference = abs(date2.tt - date1.tt)
intervals = max(1, int(days_difference / 365.25 / 5))  # 5-year intervals

# Ensure that the number of intervals is a whole number plus one for the end date
intervals = intervals + 1 if days_difference % (365.25 * 5) != 0 else intervals

dates = [date1 + (date2 - date1) * (i / intervals) for i in range(intervals + 1)]

try:
    neptune_distances = [calculate_celestial_distance(ephemeris, 'neptune barycenter', date) for date in dates]
    pluto_distances = [calculate_celestial_distance(ephemeris, 'pluto barycenter', date) for date in dates]
except Exception as e:
    print(f"Error calculating distances: {e}")
else:
    # Calculate velocities
    neptune_velocities = calculate_velocities(dates, neptune_distances)
    pluto_velocities = calculate_velocities(dates, pluto_distances)

  # Calculate the foces and accelerations
    mass_pluto = 1.30900e22
    mass_neptune = 1.024e26
    mass_sun = 1.989e30

    neptune_forces = [calculate_gravitational_force(mass_neptune, mass_sun, distance) if distance is not None else None for distance in neptune_distances]
    neptune_accelerations = [calculate_acceleration(force, mass_neptune) if force is not None else None for force in neptune_forces]

    pluto_forces = [calculate_gravitational_force(mass_pluto, mass_sun, distance) if distance is not None else None for distance in pluto_distances]
    pluto_accelerations = [calculate_acceleration(force, mass_pluto) if force is not None else None for force in pluto_forces]

    # Convert dates to years for plotting, using the first year as a reference
    year_labels = np.linspace(start=date1.utc.year, stop=date2.utc.year, num=len(dates))

    # Start of Output #
    print("\n_______________________________________")
    print("Start of Output")
    print("Rate of change is calculated by central difference method w/neighbor positions.")
    print("Acceleration is calculated by Newton's law of universal gravitation w/masses of bodies & Sun.")
    print("_______________________________________\n")

    # Print Tableable
    print_distance_table(year_labels, neptune_distances, neptune_velocities, neptune_accelerations, pluto_distances, pluto_velocities, pluto_accelerations)

    # Prints intervals on which Pluto is closer to the Sun than Neptune
    print_closer_intervals(year_labels, neptune_distances, pluto_distances)

    #Print Graph
    plt.figure(figsize=(14, 6))
    plt.plot(year_labels, neptune_distances, label='Neptune')
    plt.plot(year_labels, pluto_distances, label='Pluto')
    plt.xlabel('Year')
    plt.ylabel('Distance (AU)')
    plt.title('Distance from the Sun to Neptune and Pluto over Time')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
