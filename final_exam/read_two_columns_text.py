"""
This module can read data from a two-column text file and return the data as a numpy ndarray.
__author__ = "Jacob, Noah"
"""

import numpy as np

def read_two_columns_text(filename):
    """
    Reads the first two columns from a specified text file and returns the data as a transposed numpy ndarray.

    Parameters:
    filename (str): The path to the text file to be read.

    Returns:
    ndarray: A numpy ndarray with shape (2, M), where M is the number of data points.

    Raises:
    OSError: If the file cannot be found or read.
    """
    try:
        data = np.loadtxt(filename, usecols=(0, 1))
        return data.T
    except OSError:
        raise OSError(f"File {filename} cannot be found for reading.")

if __name__ == "__main__":
    # Example usage: Using a text file with two columns and nine rows
    filename = 'final_exam/volumes_energies.dat'
    data = read_two_columns_text(filename)
    print(f"{data=}, shape={data.shape}")
