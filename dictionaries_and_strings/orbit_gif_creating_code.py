import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Set the path for ImageMagick if it is not in the PATH environment variable
plt.rcParams['animation.convert_path'] = 'path_to_imagemagick'

# Initialization function
def init():
    trajectory.set_data([], [])
    planet.set_data([], [])
    velocity_text.set_text('')
    return trajectory, planet, velocity_text

# Animation function
def animate(i):
    trajectory.set_data(x[:i], y[:i])
    planet.set_data(x[i], y[i])
    velocity_text.set_text('Velocity magnitude = {:.2f}'.format(velocity_magnitude[i]))
    return trajectory, planet, velocity_text


# Example of showing the animation
if __name__ == '__main__':
    # Constants
    a = 1 # Semi-major axis, test value
    b = 1  # Semi-minor axis, test value for circle
    omega = 1  # Angular velocity, test value

    # Time discretization
    n_frames = 100  # Number of frames in the animation
    t = np.linspace(0, 2 * np.pi / omega, n_frames)

    # Position calculations
    x = a * np.cos(omega * t)
    y = b * np.sin(omega * t)

    # Velocity calculations
    v_x = -omega * a * np.sin(omega * t)
    v_y = omega * b * np.cos(omega * t)
    velocity_magnitude = omega * np.sqrt(a**2 * np.sin(omega * t)**2 + b**2 * np.cos(omega * t)**2)

    # Set up the figure, the axis, and the plot element to animate
    fig, ax = plt.subplots()
    ax.set_aspect('equal')  # This line ensures that the aspect ratio is equal
    planet, = ax.plot([], [], 'bo', ms=10)  # Planet as a blue circle
    trajectory, = ax.plot([], [], 'r-')  # Planet trajectory as a red line
    velocity_text = ax.text(0.02, 0.85, '', transform=ax.transAxes)
    ax.set_xlim(-a * 2.0, a * 2.0)
    ax.set_ylim(-b * 2.0, b * 2.0)

    # Call the animator
    ani = FuncAnimation(fig, animate, frames=n_frames, init_func=init, blit=True)

    # Save the animation
    ani.save('orbit_animation.gif', writer='imagemagick', fps=30)

    # Show the plot (optional, you might want to comment this out when running the script to generate the GIF)
    plt.show()
