import numpy as np

def gaussian(position):
    # Gaussian function calculation
    return (1/np.sqrt(2*np.pi)) * np.exp(-0.5 * position**2)

if __name__ == '__main__':
    # Generating 41 uniformly spaced coordinates between -4 and 4
    positions = np.linspace(-4, 4, 41)
    # Calculating Gaussian values for each position
    gaussian_values = [gaussian(x) for x in positions]
    # Creating a list of tuples with positions and their corresponding Gaussian values
    positions_and_values = list(zip(positions, gaussian_values))
    # Printing the positions and Gaussian values
    for position, value in positions_and_values:
        print(f"Position: {position:.1f}, Gaussian Value: {value:.6f}")
