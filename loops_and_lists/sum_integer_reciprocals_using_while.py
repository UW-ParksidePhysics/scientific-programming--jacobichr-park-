import math

summation = 0
starting_index = 1
index = starting_index
maximum_index = 100

while index <= maximum_index:
    summation += 1/index
    index += 1

print(f'sum(k = {starting_index}, {maximum_index}) 1/k = {summation}')

# Error 1: 'while index < maximum_index' changed to 'while index <= maximum_index' #
# Error 2: added 'index += 1' to while loop #
# Error 3: No 3rd error found, code output agrees with personal calculator #

#summation = 0
#starting_index = 1
#index = starting_index
#maximum_index = 100

#while index < maximum_index:
#    summation += 1/index

#print(f'sum(k = {starting_index}, {maximum_index}) 1/k = {summation}')