import math

initial_height = 10.0 #meters
initial_velocity = 1.0 #m/s
gravity_phobos = 0.0057 #m/s^2
gravity_deimos = 0.003 #m/s^2

time_list = []
position_list_deimos = []
position_list_phobos = []

def freefall_height(time, initial_velocity, gravity):
  return(initial_velocity * time - 0.5 * gravity * time**2)


# Header #
print("For initial velocity of 1.00 m/s:")
print("-----------------------------------------------------")
print("Phobos (g = 0.0057 m/s/s)  Deimos (g = 0.003 m/s/s)")
print("-----------------------------------------------------")
print("times_positions data")
print("Phobos\t\t Deimos")
print("t(s)         y(m)")

for time in range(0, int((2 * initial_velocity / gravity_deimos)) + 50, 50):
  time_list.append(time)
  position_list_deimos.append(freefall_height(time, initial_velocity, gravity_deimos))
  position_list_phobos.append(freefall_height(time, initial_velocity, gravity_phobos))


times_positions = list(zip(time_list, position_list_deimos))


for time, position_deimos in times_positions:
  print(f'{time:.2f} (s) \t {position_deimos:.2f} (m)')
  

print("-----------------------------------------------------")
print("time_positions data")
print("Phobos\t\t Deimos")
print("t(s)         y(m)")

time_positions = [
  (0.00, 0.00),
  (50.00, 46.25),
  (100.00, 85.00),
  (150.00, 116.25),
  (200.00, 140.00),
  (250.00, 156.25),
  (300.00, 165.00),
  (350.00, 166.25),
  (400.00, 160.00),
  (450.00, 146.25),
  (500.00, 125.00),
  (550.00, 96.25),
  (600.00, 60.00),
  (650.00, 16.25),
  (700.00, -35.00)
]
# Loop over the time_positions list and write out the t and y values with two decimal places
for t, y in time_positions:
  print(f'{t:.2f} (s) \t {y:.2f} (m)')