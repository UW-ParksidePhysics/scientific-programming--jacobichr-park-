# Initialize the nearby star data list
nearby_star_data = [
    ('Alpha Centauri A',    4.3,  0.26,      1.56),
    ('Alpha Centauri B',    4.3,  0.077,     0.45),
    ('Alpha Centauri C',    4.2,  0.00001,   0.00006),
    ("Barnard's Star",      6.0,  0.00004,   0.0005),
    ('Wolf 359',            7.7,  0.000001,  0.00002),
    ('BD +36 degrees 2147', 8.2,  0.0003,    0.006),
    ('Luyten 726-8 A',      8.4,  0.000003,  0.00006),
    ('Luyten 726-8 B',      8.4,  0.000002,  0.00004),
    ('Sirius A',            8.6,  1.00,      23.6),
    ('Sirius B',            8.6,  0.001,     0.003),
    ('Ross 154',            9.4,  0.00002,   0.0005), 
]

def print_table(data, header):
    print(f"{header[0]:<25} {header[1]}")
    for item in data:
        print(f"{item[0]:<25} {item[1]}")

# Sort the data by distance
sorted_by_distance = sorted(nearby_star_data, key=lambda x: x[1])
print("Star Name vs. Distance:")
print_table(sorted_by_distance, ("Star Name", "Distance (ly)"))

# Sort the data by apparent brightness
sorted_by_apparent_brightness = sorted(nearby_star_data, key=lambda x: x[2], reverse=True)
print("\nStar Name vs. Apparent Brightness:")
print_table(sorted_by_apparent_brightness, ("Star Name", "Apparent Brightness"))

# Sort the data by luminosity
sorted_by_luminosity = sorted(nearby_star_data, key=lambda x: x[3], reverse=True)
print("\nStar Name vs. Luminosity:")
print_table(sorted_by_luminosity, ("Star Name", "Luminosity"))