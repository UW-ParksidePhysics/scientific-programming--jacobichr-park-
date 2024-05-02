"""
This module fits a quadratic polynomial to an ndarray of shape (2, M) using NumPy's polyfit function.
__author__ = "Jacob, Noah"
"""

import numpy as np

def calculate_quadratic_fit(data, threshold=1e-10):
    """
    Fits a quadratic polynomial to provided x-y data and returns the coefficients,
    rounding small coefficients (below a given threshold) to zero.

    Parameters:
    data (ndarray): A numpy ndarray with shape (2, M), where 2 represents x-y data and M is the number of data points.
    threshold (float): The threshold below which coefficients are rounded to zero.

    Returns:
    ndarray: An array of shape (3) containing the polynomial coefficients: [constant, linear, quadratic].

    Raises:
    IndexError: If the data array has inappropriate dimensions.
    """
    if data.shape[0] != 2 or data.shape[1] < 3:
        raise IndexError("Data array must have 2 rows and at least 3 data points to fit a quadratic polynomial.")

    # Fit a second-degree polynomial to the data
    coefficients = np.polyfit(data[0], data[1], 2)
    # Reverse the order to match the requested output: [constant, linear, quadratic]
    reversed_coefficients = coefficients[::-1]
    # Round near zero values based on threshold
    rounded_coefficients = np.where(np.abs(reversed_coefficients) < threshold, 0, reversed_coefficients)
    return rounded_coefficients

if __name__ == "__main__":
    # Generate test data: y = x^2 (should fit exactly to a quadratic equation with coefficients [0, 0, 1])
    x = np.linspace(-1, 1, 100)
    y = x**2
    data = np.vstack((x, y))

    # Calculate quadratic fit:
    quadratic_coefficients = calculate_quadratic_fit(data)
    print("Quadratic Coefficients:")
    print(quadratic_coefficients)
