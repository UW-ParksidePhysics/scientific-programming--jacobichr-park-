import math

float_list = []
starting_value = 1.0
final_value = 2.0
step = 0.05

for current_value in [value / 100 for value in range(int(starting_value * 100), int(final_value * 100) + 1, int(step * 100))]:
    float_list.append(current_value)


float_list_comprehension = [value / 100 for value in range(int(starting_value * 100), int(final_value * 100) + 1, int(step * 100))]

print(f'Float list for loop: {float_list}')
print(f'Float list comprehension: {float_list_comprehension}')
