import vpython as vp

# Initialize parameters for the first ball
initial_position1 = vp.vector(-5., 0., 1.)
initial_velocity1 = vp.vector(1., 0.2, -0.2)

# Create the first ball
ball1 = vp.sphere(pos=initial_position1, radius=0.1, color=vp.color.red, make_trail=True)

# Initialize parameters for the second ball
initial_position2 = vp.vector(0., -5., 1.)
initial_velocity2 = vp.vector(0.2, 1., -0.2)

# Create the second ball
ball2 = vp.sphere(pos=initial_position2, radius=0.1, color=vp.color.blue, make_trail=True)

# Set up the animation parameters
animation_time_step = 0.1  # seconds
rate_of_animation = 1 / animation_time_step
time_step = 0.05
stop_time = 10.

time = 0.
while time < stop_time:
    vp.rate(rate_of_animation)
    # Update positions of ball1
    x1 = initial_position1.x + initial_velocity1.x * time
    y1 = initial_position1.y + initial_velocity1.y * time
    z1 = initial_position1.z + initial_velocity1.z * time
    ball1.pos = vp.vector(x1, y1, z1)

    # Update positions of ball2
    x2 = initial_position2.x + initial_velocity2.x * time
    y2 = initial_position2.y + initial_velocity2.y * time
    z2 = initial_position2.z + initial_velocity2.z * time
    ball2.pos = vp.vector(x2, y2, z2)

    time += time_step
