from vpython import *
import numpy as np

# Gravitational accelerations (m/s^2)
g_earth = 9.81
g_mars = 3.71
g_moon = 1.62

# Launch angles for two different scenarios
theta1 = np.pi / 4  # 45 degrees

# Function to calculate initial velocity components based on the given equations
def initial_velocity_components(g, theta):
    # Assuming we want the projectile to land at the opposite corner of a 100m field
    R = 100  # Range R (meters)
    # Time of flight equation T = 2*v0y / ay, solved for v0y to achieve the desired range
    v0y = sqrt(R * g / sin(2 * theta))
    v0x = v0y / tan(theta)  # Horizontal velocity component
    return v0x, v0y

# Function to create and launch the ball
def launch_ball(field_pos, color, g, theta):
    # Initial velocity
    v0x, v0y = initial_velocity_components(g, theta)

    # Create ball with trail
    ball = sphere(pos=field_pos + vec(-50, 0, 0), radius=1, color=color, make_trail=True, trail_radius=0.5)
    ball.velocity = vec(v0x, v0y, 0)

    # Time of flight based on the initial vertical velocity and gravitational acceleration
    T = 2 * v0y / g

    return ball, T

# Set up the display
scene = canvas(title='Projectile Motion Comparison')

# Create fields
field_earth = box(pos=vec(-150, -5, 0), size=vec(100, 0.1, 50), color=color.green)
field_mars = box(pos=vec(0, -5, 0), size=vec(100, 0.1, 50), color=color.red)
field_moon = box(pos=vec(150, -5, 0), size=vec(100, 0.1, 50), color=color.gray)

# Launch balls for each celestial body and get their time of flight
ball_earth, T_earth = launch_ball(field_earth.pos, color.blue, g_earth, theta1)
ball_mars, T_mars = launch_ball(field_mars.pos, color.pink, g_mars, theta1)
ball_moon, T_moon = launch_ball(field_moon.pos, color.white, g_moon, theta1)

# Maximum time of flight for stopping the simulation
T_max = max(T_earth, T_mars, T_moon)

# Time increment
dt = min(T_earth, T_mars, T_moon) / 300
t = 0

# Simulation loop
while t < T_max:
    rate(300)  # Limit the number of loops per second

    # Update positions and velocities
    if ball_earth.pos.y >= 0:
        ball_earth.pos += ball_earth.velocity * dt
        ball_earth.velocity.y -= g_earth * dt

    if ball_mars.pos.y >= 0:
        ball_mars.pos += ball_mars.velocity * dt
        ball_mars.velocity.y -= g_mars * dt

    if ball_moon.pos.y >= 0:
        ball_moon.pos += ball_moon.velocity * dt
        ball_moon.velocity.y -= g_moon * dt

    # Update time
    t += dt

# Add text labels for each celestial body's field
label(pos=field_earth.pos - vec(0, 10, 0), text='EARTH', color=color.green, height=10, box=False)
label(pos=field_mars.pos - vec(0, 10, 0), text='MARS', color=color.red, height=10, box=False)
label(pos=field_moon.pos - vec(0, 10, 0), text='MOON', color=color.gray, height=10, box=False)
