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

def print_table(data, headers):
  format_str = "{:<25} {:>15}"  # Adjusted format string for two columns
  print(format_str.format(*headers))
  for item in data:
      print(format_str.format(item[0], item[1]))

if __name__ == "__main__":
  # Sort the data by distance and print the table
  sorted_by_distance = sorted(nearby_star_data, key=lambda x: x[1], reverse=True)
  print("Star Name vs. Distance:")
  print_table([(star[0], star[1]) for star in sorted_by_distance], ("Star Name", "Distance (ly)"))

  # Sort the data by apparent brightness and print the table
  sorted_by_apparent_brightness = sorted(nearby_star_data, key=lambda x: x[2], reverse=True)
  print("\nStar Name vs. Apparent Brightness:")
  print_table([(star[0], star[2]) for star in sorted_by_apparent_brightness], ("Star Name", "Apparent Brightness"))

  # Sort the data by luminosity and print the table
  sorted_by_luminosity = sorted(nearby_star_data, key=lambda x: x[3], reverse=True)
  print("\nStar Name vs. Luminosity:")
  print_table([(star[0], star[3]) for star in sorted_by_luminosity], ("Star Name", "Luminosity"))
