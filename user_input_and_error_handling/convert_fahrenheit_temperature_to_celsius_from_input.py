# Python program to convert temperature from Fahrenheit to Celsius

def fahrenheit_to_celsius(fahrenheit):
    return (fahrenheit - 32) * 5.0/9.0

def print_tempurature_conversion():
    # Ask the user for the temperature in Fahrenheit
    fahrenheit = float(input("Enter temperature in Fahrenheit: "))

    # Convert the temperature to Celsius
    celsius = fahrenheit_to_celsius(fahrenheit)

    # Print the temperatures in both Fahrenheit and Celsius
    print(f"{fahrenheit:.2f} Fahrenheit is equivalent to {celsius:.2f} Celsius.")

if __name__ == "__main__":
  print_tempurature_conversion()
