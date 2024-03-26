import sys

# Python program to convert temperature from Fahrenheit to Celsius

def fahrenheit_to_celsius(fahrenheit):
    return (fahrenheit - 32) * 5.0 / 9.0

def print_temperature_conversion():
    try:
        # Check if the temperature was passed as a command line argument
        if len(sys.argv) != 2:
            raise ValueError("Temperature in Fahrenheit is required as a command-line argument.")

        # Read the temperature in Fahrenheit from the command line argument
        fahrenheit = float(sys.argv[1])

        # Convert the temperature to Celsius
        celsius = fahrenheit_to_celsius(fahrenheit)

        # Print the temperatures in both Fahrenheit and Celsius
        print(f"{fahrenheit:.2f} Fahrenheit is equivalent to {celsius:.2f} Celsius.")
    except ValueError as e:
        print(f"Error: {e}")
        print("Usage: python3 convert_temperature_from_command_line.py <temperature_in_Fahrenheit>")
        print("You must input as a command-line argument, as either an integer or a float!")
        sys.exit(1)

if __name__ == "__main__":
    print_temperature_conversion()
