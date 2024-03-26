def calculate_velocity_and_acceleration(positions, times):
    if len(positions) < 2:
        raise ValueError("Need at least two positions to calculate velocity.")

    velocities = []
    accelerations = [0]  # Acceleration is zero for constant velocity

    for i in range(1, len(positions)):
        velocity = (positions[i] - positions[i - 1]) / (times[i] - times[i - 1])
        velocities.append(velocity)

    return velocities[0], accelerations[0]  # Return the constant velocity and zero acceleration

def test_kinematics():
    # Define the times and initial velocity
    times = [0, 0.5, 1.5, 2.2]
    v = 2  # Assuming a constant velocity of 2 units/time

    # Calculate positions based on constant velocity
    positions = [v * t for t in times]

    print("Testing kinematics for constant velocity...")

    velocity, acceleration = calculate_velocity_and_acceleration(positions, times)
    for i in range(len(times)):
        print(f"\nTime: {times[i]}")
        print(f"Input Position: {positions[i]}")
        print(f"Calculated Velocity: {velocity}")
        print(f"Calculated Acceleration: {acceleration}")

if __name__ == "__main__":
  test_kinematics()
