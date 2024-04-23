import vpython as vp

# Define the balls with different initial positions and velocities
ball1 = vp.sphere(pos=vp.vector(-10., -7., 0.), radius=0.5, color=vp.color.blue, make_trail=True)
ball2 = vp.sphere(pos=vp.vector(-10., 7., 0.), radius=0.5, color=vp.color.green, make_trail=True)

initial_velocity1 = vp.vector(25., 15., 0.)  # Ball 1 (blue) moving at an angle towards the wall
initial_velocity2 = vp.vector(25., -5., 0.)  # Ball 2 (green) moving at an angle towards the wall

# Define the wall
wall_center = vp.vector(0., 0., 0.)
wall_dimensions = vp.vector(0.25, 10., 10.)
wall = vp.box(pos=wall_center, size=wall_dimensions, color=vp.color.red)

# Animation parameters
animation_time_step = 0.01  # seconds
rate_of_animation = 1 / animation_time_step
time_step = 0.005
stop_time = 1.

time = 0.
while time < stop_time:
    vp.rate(rate_of_animation)

    # Update positions of ball 1
    ball1.pos.x += initial_velocity1.x * time_step
    ball1.pos.y += initial_velocity1.y * time_step
    # Check collision for ball 1
    if abs(ball1.pos.x - wall.pos.x) <= (ball1.radius + wall.size.x / 2):
        initial_velocity1.x = -initial_velocity1.x  # reverse ball 1 velocity

    # Update positions of ball 2
    ball2.pos.x += initial_velocity2.x * time_step
    ball2.pos.y += initial_velocity2.y * time_step
    # Check collision for ball 2
    if abs(ball2.pos.x - wall.pos.x) <= (ball2.radius + wall.size.x / 2):
        initial_velocity2.x = -initial_velocity2.x  # reverse ball 2 velocity

    time += time_step
