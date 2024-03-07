import math

maximum_integer = 10
for_loop_number = 0

for i in range(1, maximum_integer + 1):
  for_loop_number = for_loop_number + i

print(f'The sum from a for loop:{for_loop_number}')
print(f'The equation sum:{maximum_integer * (maximum_integer + 1) / 2}')