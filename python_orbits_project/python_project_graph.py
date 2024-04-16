from skyfield.api import load, Loader
import matplotlib.pyplot as plt
import numpy as np

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
    # Convert dates to years for plotting, using the first year as a reference
    year_labels = np.linspace(start=date1.utc.year, stop=date2.utc.year, num=len(dates))

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
