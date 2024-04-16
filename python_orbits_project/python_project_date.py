from skyfield.api import load, Topos
from datetime import datetime

def is_valid_date(day, month, year):
    try:
        datetime(year, month, day)
        return True
    except ValueError:
        return False

def calculate_celestial_distance(celestial_body_name, day, month, year):
    planets = load('de421.bsp')
    sun = planets['sun']
    celestial_body = planets[celestial_body_name]

    ts = load.timescale()
    t = ts.utc(year, month, day)

    distance = sun.at(t).observe(celestial_body).distance().au  # Removed .apparent() for direct distance

    return distance

# Get date input from user
day = int(input("Enter the day (DD): "))
month = int(input("Enter the month (MM): "))
year = int(input("Enter the year (YYYY): "))

if is_valid_date(day, month, year):
    neptune_distance = calculate_celestial_distance('neptune barycenter', day, month, year)
    pluto_distance = calculate_celestial_distance('pluto barycenter', day, month, year)

    print(f"On {month:02d}/{day:02d}/{year}, the distance of Neptune from the Sun is approximately {neptune_distance:.2f} AU.")
    print(f"On {month:02d}/{day:02d}/{year}, the distance of Pluto from the Sun is approximately {pluto_distance:.2f} AU.")
else:
    print("The entered date is invalid.")
