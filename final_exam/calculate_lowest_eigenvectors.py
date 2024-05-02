"""
This module calculates the eigenvectors with the smallest eigenvalues for a given square matrix using NumPy's eig function.
__author__ = "Jacob, Noah"
"""

import numpy as np

def calculate_lowest_eigenvectors(square_matrix, number_of_eigenvectors=3):
    """
    Finds the eigenvectors associated with the lowest eigenvalues of a square matrix.

    Parameters:
    square_matrix (ndarray): A square matrix of shape (M, M).
    number_of_eigenvectors (int): The number of eigenvectors to return. Default is 3.

    Returns:
    eigenvalues (ndarray): An array of the K lowest eigenvalues, sorted from lowest to highest.
    eigenvectors (ndarray): An array of shape (K, M) containing the K eigenvectors corresponding to the smallest eigenvalues.
    """
    if square_matrix.shape[0] != square_matrix.shape[1]:
        raise ValueError("Input must be a square matrix.")
    if number_of_eigenvectors > square_matrix.shape[0]:
        raise ValueError("Number of eigenvectors requested exceeds matrix dimensions.")

    # Calculate eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(square_matrix)

    # Sort eigenvalues and select the smallest ones
    idx = eigenvalues.argsort()[:number_of_eigenvectors]
    sorted_eigenvalues = eigenvalues[idx]
    sorted_eigenvectors = eigenvectors[:, idx].T

    return sorted_eigenvalues, sorted_eigenvectors

if __name__ == "__main__":
    # Test the function with a specified matrix and number of eigenvectors
    square_matrix = np.array([[2, -1], [-1, 2]])
    number_of_eigenvectors = 2

    eigenvalues, eigenvectors = calculate_lowest_eigenvectors(square_matrix, number_of_eigenvectors)
    print("Eigenvalues:")
    print(eigenvalues)
    print("Eigenvectors:")
    print(eigenvectors)
