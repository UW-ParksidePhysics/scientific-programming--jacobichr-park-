import math

fahrenheit = 10
while fahrenheit < 100:
  celsius = (fahrenheit - 32) * 5/9
  approximate = (fahrenheit-30)/2
  print(f"{fahrenheit}°F, |  {celsius:.2f}°C, | Approximatly {approximate:.2f}°C ")
  fahrenheit += 10