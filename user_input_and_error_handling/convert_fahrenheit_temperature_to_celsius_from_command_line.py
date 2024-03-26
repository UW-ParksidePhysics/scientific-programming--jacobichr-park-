import sys

# Python program to convert temperature from Fahrenheit to Celsius

def fahrenheit_to_celsius(fahrenheit):
    return (fahrenheit - 32) * 5.0/9.0

def print_tempurature_conversion():
    # Check if the temperature was passed as a command line argument
    if len(sys.argv) != 2:
        print("Usage: python program.py <temperature_in_Fahrenheit>")
        sys.exit(1)

    # Read the temperature in Fahrenheit from the command line argument
    fahrenheit = float(sys.argv[1])

    # Convert the temperature to Celsius
    celsius = fahrenheit_to_celsius(fahrenheit)

    # Print the temperatures in both Fahrenheit and Celsius
    print(f"{fahrenheit:.2f} Fahrenheit is equivalent to {celsius:.2f} Celsius.")

if __name__ == "__main__":
  print_tempurature_conversion()
