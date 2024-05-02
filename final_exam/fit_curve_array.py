"""
This module generates x-y data points based on quadratic coefficients over a specified range and number of points.
Outputs points in scientific (e+...) notation.
__author__ = "Jacob, Noah"
"""

import numpy as np

def fit_curve_array(quadratic_coefficients, minimum_x, maximum_x, number_of_points=100):
    """
    Generates a curve from quadratic coefficients over a specified x range with a given number of points.

    Parameters:
    quadratic_coefficients (ndarray): An array of shape (3) containing the quadratic coefficients [constant, linear, quadratic].
    minimum_x (float): The minimum x-value for the curve.
    maximum_x (float): The maximum x-value for the curve.
    number_of_points (int): The number of points N to evaluate the quadratic function. Default is 100.

    Returns:
    ndarray: An array of shape (2, N) where 2 represents x-y data and N is the number of evaluation points.

    Raises:
    ArithmeticError: If maximum_x < minimum_x.
    IndexError: If number_of_points <= 2.
    """
    if maximum_x < minimum_x:
        raise ArithmeticError("maximum_x must be greater than minimum_x.")
    if number_of_points <= 2:
        raise IndexError("number_of_points must be greater than 2.")

    # Generate x values
    x_values = np.linspace(minimum_x, maximum_x, number_of_points)
    # Evaluate the polynomial
    y_values = np.polyval(quadratic_coefficients[::-1], x_values)  # Ensure coefficients are in [quadratic, linear, constant] order for polyval

    # Stack x and y values to create the fit curve
    fit_curve = np.vstack((x_values, y_values))
    return fit_curve

if __name__ == "__main__":
    # Test the module with a quadratic equation y = x^2
    quadratic_coefficients = [0, 0, 1]  # [constant, linear, quadratic]
    minimum_x = -2
    maximum_x = 2
    fit_curve = fit_curve_array(quadratic_coefficients, minimum_x, maximum_x)
    print("Fit Curve:")
    print(fit_curve)
    print("Shape of fit curve array:", fit_curve.shape)
