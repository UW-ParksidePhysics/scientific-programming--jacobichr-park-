import math

def convert_celsius_to_fahrenheit(celsius_temperature): # Converts input temperature from Celsius to Fahrenheit
    return (9./5) * celsius_temperature + 32


def convert_fahrenheit_to_celsius(fahrenheit_temperature): # Converts input temperature from Fahrenheit to Celsius
  return (5./9) * (fahrenheit_temperature - 32)


if __name__ == "__main__":
  celsius_tempuratures = [0,21,100]
  for some_celsius_temperature in celsius_tempuratures:
    converted_tempurature = convert_celsius_to_fahrenheit(some_celsius_temperature)
    print(f"{some_celsius_temperature:.2f} degrees Celsius is {converted_tempurature:.2f} degrees Fahrenheit which is stil {convert_fahrenheit_to_celsius(converted_tempurature):.2f} degrees Celsius!")