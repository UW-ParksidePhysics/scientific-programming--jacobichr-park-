# Input:
# import user_input_and_error_handling.convert_temperature as convert_temperature

# print(convert_temperature.test_conversion())
# print(convert_temperature.handle_conversion(21.3, 'C'))
# print(convert_temperature.handle_conversion(21.3, 'F'))
# print(convert_temperature.kelvin_to_celsius(100))

# Output:
# True
# 21.30 C is 70.34 F or 294.45 K
# None
# 21.30 F is -5.94 C or 267.21 K
# None
# -173.14999999999998


import sys

def celsius_to_fahrenheit(celsius):
    return celsius * 9.0 / 5.0 + 32

def fahrenheit_to_celsius(fahrenheit):
    return (fahrenheit - 32) * 5.0 / 9.0

def celsius_to_kelvin(celsius):
    return celsius + 273.15

def kelvin_to_celsius(kelvin):
    return kelvin - 273.15

def fahrenheit_to_kelvin(fahrenheit):
    return (fahrenheit - 32) * 5.0 / 9.0 + 273.15

def kelvin_to_fahrenheit(kelvin):
    return (kelvin - 273.15) * 9.0 / 5.0 + 32

def test_conversion():
  success = True
  tolerance = 0.01

  # Sample temperatures for testing
  fahrenheit_temperature = 32.0  # Freezing point of water in Fahrenheit
  celsius_temperature = 0.0      # Freezing point of water in Celsius
  kelvin_temperature = 273.15    # Freezing point of water in Kelvin

  # Test Fahrenheit to Celsius and back
  if abs(fahrenheit_to_celsius(celsius_to_fahrenheit(fahrenheit_temperature)) - fahrenheit_temperature) > tolerance:
      print("Fahrenheit to Celsius back to Fahrenheit conversion failed.")
      success = False

  # Test Celsius to Kelvin and back
  if abs(celsius_to_kelvin(kelvin_to_celsius(celsius_temperature)) - celsius_temperature) > tolerance:
      print("Celsius to Kelvin back to Celsius conversion failed.")
      success = False

  # Test Kelvin to Fahrenheit and back
  if abs(kelvin_to_fahrenheit(fahrenheit_to_kelvin(kelvin_temperature)) - kelvin_temperature) > tolerance:
      print("Kelvin to Fahrenheit back to Kelvin conversion failed.")
      success = False

  # Using assert to indicate a failure in conversion tests
  assert success, "Not all temperature conversion tests passed."

  return success


def handle_conversion(temperature, scale):
    if scale.upper() == 'C':
        fahrenheit = celsius_to_fahrenheit(temperature)
        kelvin = celsius_to_kelvin(temperature)
        print(f"{temperature:.2f} C is {fahrenheit:.2f} F or {kelvin:.2f} K")
    elif scale.upper() == 'F':
        celsius = fahrenheit_to_celsius(temperature)
        kelvin = fahrenheit_to_kelvin(temperature)
        print(f"{temperature:.2f} F is {celsius:.2f} C or {kelvin:.2f} K")
    elif scale.upper() == 'K':
        celsius = kelvin_to_celsius(temperature)
        fahrenheit = kelvin_to_fahrenheit(temperature)
        print(f"{temperature:.2f} K is {celsius:.2f} C or {fahrenheit:.2f} F")
    else:
        print("Invalid scale. Please use C, F, or K.")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == "verify":
            test_conversion()
            print("All conversion tests passed successfully.")
        elif len(sys.argv) == 3:
            try:
                temp = float(sys.argv[1])
                scale = sys.argv[2]
                handle_conversion(temp, scale)
            except ValueError:
                print("Invalid temperature. Please enter a valid number.")
        else:
            print("Invalid usage.")
            print("Usage: program.py <temperature> <scale>")
            print("Example: program.py 21.3 C")
            print("Or use 'verify' to test the conversion functions.")
    else:
        print("Usage: program.py <temperature> <scale>")
        print("Example: program.py 21.3 C")
        print("Or use 'verify' to test the conversion functions.")
