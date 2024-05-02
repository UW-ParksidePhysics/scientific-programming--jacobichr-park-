"""
This module can read data from a Numpy ndarray and preform an analysis of
various statistical properties of the data.
__author__ = "Jacob, Noah"
"""

import numpy as np
from scipy import stats

def calculate_bivariate_statistics(data):
    """
    Calculates statistical characteristics for bivariate (x, y) data using scipy's stats.describe function.

    Parameters:
    data (ndarray): A numpy ndarray with shape (2, M), where 2 represents x-y data and M is the number of data points.

    Returns:
    ndarray: An array of shape (6) containing statistics:
        - mean of y (statistics[0])
        - standard deviation of y (statistics[1])
        - minimum x-value (statistics[2])
        - maximum x-value (statistics[3])
        - minimum y-value (statistics[4])
        - maximum y-value (statistics[5])

    Raises:
    IndexError: If the data array has inappropriate dimensions.
    """
    if data.shape[0] != 2 or data.shape[1] <= 1:
        raise IndexError("Data array must have 2 rows and more than 1 column.")

    stats_y = stats.describe(data[1, :])
    min_x = np.min(data[0, :])
    max_x = np.max(data[0, :])
    min_y = stats_y.minmax[0]
    max_y = stats_y.minmax[1]

    return np.array([
        stats_y.mean, np.sqrt(stats_y.variance),
        min_x, max_x, min_y, max_y
    ])

if __name__ == "__main__":
    # Example usage: Using an array centered around zero of (x ∈ [-10, 10])
    # Create array:
    x = np.linspace(-10, 10, 101)
    y = x**2
    data = np.vstack((x, y))

    # Preform analysis:
    statistics = calculate_bivariate_statistics(data)
    print("Statistics:")
    print("Mean of Y:", statistics[0])
    print("Standard Deviation of Y:", statistics[1])
    print("Minimum X-value:", statistics[2])
    print("Maximum X-value:", statistics[3])
    print("Minimum Y-value:", statistics[4])
    print("Maximum Y-value:", statistics[5])
