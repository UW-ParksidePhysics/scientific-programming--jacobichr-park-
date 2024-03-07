import math

initial_height = 10.0 #meters
initial_velocity = 1.0 #m/s
gravity_phobos = 0.0057 #m/s^2
gravity_deimos = 0.003 #m/s^2

def freefall_height(time, initial_velocity, gravity):
  return(initial_velocity * time - 0.5 * gravity * time**2)


# Header #
print("For initial velocity of 1.00 m/s:")
print("-----------------------------------------------------")
print("Phobos (g = 0.0057 m/s/s)  Deimos (g = 0.003 m/s/s)")
print("-----------------------------------------------------")
print("Using a for loop:")
print("Phobos\t\t\t Deimos")
print("t(s)     y(m)    t(s)    y(m)")

for time in range(0, int((2 * initial_velocity / gravity_deimos)) + 50, 50):
  if freefall_height(time, initial_velocity, gravity_phobos) >= 0:
    print(f'{time:<5} \t {freefall_height(time, initial_velocity, gravity_phobos):.2f} \t {time:<5} \t {freefall_height(time, initial_velocity, gravity_deimos):.2f}')
  else:
    print(f'\t\t\t\t {time:<5} \t{freefall_height(time, initial_velocity, gravity_deimos):.2f}')

print("-----------------------------------------------------")
print("Using a while loop:")
print("Phobos\t\t\t Deimos")
print("t(s)     y(m)    t(s)    y(m)")
time = 0
while time < int((2 * initial_velocity / gravity_deimos)) + 50:
    if freefall_height(time, initial_velocity, gravity_phobos) >= 0:
        print(f'{time:<5} \t {freefall_height(time, initial_velocity, gravity_phobos):.2f} \t {time:<5} \t {freefall_height(time, initial_velocity, gravity_deimos):.2f}')
    else:
        print(f'\t\t\t\t {time:<5} \t{freefall_height(time, initial_velocity, gravity_deimos):.2f}')
    time += 50