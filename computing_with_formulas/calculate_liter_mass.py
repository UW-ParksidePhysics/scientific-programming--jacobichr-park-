iron_density = 7.87
air_density = 1.225
gasoline_density = 0.755
ice_deniity = 0.02
human_density = 1.1
silver_density = 10.49
platinum_density = 21.45

volume = 1

iron_mass = iron_density * volume
air_mass = air_density * volume
gasoline_mass = gasoline_density * volume
ice_mass = ice_deniity * volume
human_mass = human_density * volume
silver_mass = silver_density * volume
platinum_mass = platinum_density * volume

print(f'Iron mass = {iron_mass}g/cm^3')
print(f'Air mass = {air_mass}g/cm^3')
print(f'Gasoline mass = {gasoline_mass}g/cm^3')
print(f'Ice mass = {ice_mass}g/cm^3')
print(f'Human mass = {human_mass}g/cm^3')
print(f'Silver mass = {silver_mass}g/cm^3')
print(f'Platinum mass = {platinum_mass}g/cm^3')
