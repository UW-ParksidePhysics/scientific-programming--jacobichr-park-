from typing import final


interest_rate = 4.25 #Feb 12th 2024
years = 3
p = interest_rate #% per year
n = years
a = 1000 #initial amount
final_amount = a(1 + (p/100))**n

print(a)
print(final_amount)
print(final_amount-a)