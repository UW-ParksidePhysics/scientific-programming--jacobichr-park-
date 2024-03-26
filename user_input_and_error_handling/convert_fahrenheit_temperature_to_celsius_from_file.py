# Python program to convert temperature from Fahrenheit to Celsius read from a file

def fahrenheit_to_celsius(fahrenheit):
    return (fahrenheit - 32) * 5.0/9.0

def read_temperature_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if "Fahrenheit degrees:" in line:
                # Extract the temperature from the line
                temperature_str = line.split(":")[1].strip()
                return float(temperature_str)
    return None

def print_tempurature_conversion():
    # Filename where the Fahrenheit temperature is stored
    filename = 'user_input_and_error_handling/temperature_data.txt'

    # Read the Fahrenheit temperature from the file
    fahrenheit = read_temperature_from_file(filename)
    if fahrenheit is None:
        print("Temperature data not found in the file.")
        return

    # Convert the temperature to Celsius
    celsius = fahrenheit_to_celsius(fahrenheit)

    # Print the temperatures in both Fahrenheit and Celsius
    print(f"{fahrenheit:.2f} Fahrenheit is equivalent to {celsius:.2f} Celsius.")

if __name__ == "__main__":
  print_tempurature_conversion()
