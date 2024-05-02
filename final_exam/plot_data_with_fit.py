"""
This module plots x-y data points along with a fitted curve (scatter and curve plot) using matplotlib's Pyplot.
__author__ = "Jacob, Noah"
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_data_with_fit(data, fit_curve, data_format='o', fit_format='-'):
    """
    Plots given data points and a fitted curve with optional formatting.

    Parameters:
    data (ndarray): An ndarray of shape (2, M) containing the x-y data points.
    fit_curve (ndarray): An ndarray of shape (2, N) containing the x-y data for the fitted curve.
    data_format (str): Formatting string for the scatter plot data points. Default is 'o'.
    fit_format (str): Formatting string for the curve of the fit function. Default is '-'.

    Returns:
    list: A list of Line2D objects representing the plotted data and fit curve.
    """
    # Plot the scatter data
    data_line = plt.plot(data[0], data[1], data_format, label='Data Points')

    # Plot the fitted curve
    fit_line = plt.plot(fit_curve[0], fit_curve[1], fit_format, label='Fit Curve')

    # Adding legend
    plt.legend()

    # Return the list of Line2D objects created by plot
    return data_line + fit_line

if __name__ == "__main__":
    # Test the function
    data = np.array([[-2, -1, 0, 1, 2], [4, 1, 0, 1, 4]])
    x_vals = np.linspace(-2, 2, 100)
    fit_curve = np.array([x_vals, x_vals**2])

    # Visualize the data and fit
    plt.figure(figsize=(8, 6))
    combined_plot = plot_data_with_fit(data, fit_curve, data_format='x', fit_format='--')
    plt.title("Data Points with Quadratic Fit")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.show()
