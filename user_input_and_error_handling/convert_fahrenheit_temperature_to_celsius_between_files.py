def fahrenheit_to_celsius(fahrenheit):
  return (fahrenheit - 32) * 5.0 / 9.0

def read_temperatures_from_file(filename):
  temperatures = []
  with open(filename, 'r') as file:
      for line in file:
          if "Fahrenheit degrees:" in line:
              parts = line.split(":")
              if len(parts) == 2:
                  temperatures.append(float(parts[1].strip()))
  return temperatures

def write_temperatures_to_file(temperatures, filename):
  with open(filename, 'w') as file:
      file.write("Fahrenheit  Celsius\n")
      file.write("------------------\n")
      for fahrenheit in temperatures:
          celsius = fahrenheit_to_celsius(fahrenheit)
          file.write(f"{fahrenheit:11.2f}  {celsius:8.2f}\n")

def print_tempurature_conversion():
  input_filename = 'user_input_and_error_handling/multiple_temperature_datas.txt'
  output_filename = 'converted_temperatures.txt'

  # Read the Fahrenheit temperatures from the file
  fahrenheit_temperatures = read_temperatures_from_file(input_filename)

  # Write the Fahrenheit and Celsius temperatures to the output file
  write_temperatures_to_file(fahrenheit_temperatures, output_filename)
  print(f"Temperatures have been converted and written to {output_filename}")

if __name__ == "__main__":
  print_tempurature_conversion()
