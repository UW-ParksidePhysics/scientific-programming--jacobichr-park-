def parse_constants_file(filepath):
  constants = {}
  with open(filepath, 'r') as file:
      for i, line in enumerate(file):
          # Skip the header lines
          if i < 2:
              continue
          # Parse the constant name, value, and dimension
          parts = line.split()
          name = ' '.join(parts[:-2])  # The name can be multiple words
          value = float(parts[-2])  # Convert the value to float
          dimension = parts[-1]  # The dimension is the last part
          constants[name] = (value, dimension)
  return constants

if __name__ == '__main__':
  filepath = 'dictionaries_and_strings/constants.txt'
  constants_dict = parse_constants_file(filepath)
  for name, (value, dimension) in constants_dict.items():
      print(f"{name:30s}: {value} [{dimension}]")
