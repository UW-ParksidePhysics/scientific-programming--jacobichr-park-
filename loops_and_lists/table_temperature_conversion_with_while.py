import math

fahrenheit = 10
while fahrenheit < 100:
  celsius = (fahrenheit - 32) * 5/9
  print(f"{fahrenheit}° Fahrenheit is {celsius}° Celsius")
  fahrenheit += 10