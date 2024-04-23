import numpy as np
import matplotlib.pyplot as plt


# Set the variable matrix_dimension to (10)
matrix_dimension = 10

# Adjust the scale factor
s = 1 / (2 * (1 / (matrix_dimension + 1))**2)

# Create the array for the diagonal values using the new matrix_dimension
values = np.array([2 - matrix_dimension**2] + [-12 - i*10 for i in range(1, matrix_dimension)])

# Construct the 10x10 matrix H
H = s * np.diagflat(values)

# Find the eigenvalues and eigenvectors of the new H matrix
eigenvalues, eigenvectors = np.linalg.eig(H)

# Set the x_values for plotting the fifth eigenvector of the 10x10 matrix
x_values = np.linspace(1 / (matrix_dimension + 1), matrix_dimension / (matrix_dimension + 1), matrix_dimension)
fifth_eigenvector = eigenvectors[:, 4]

# Define the range for the sin function plot
x_range = np.linspace(0, 1, 1000)
sin_function = np.sqrt(2) * np.sin(np.pi * x_range)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_values, fifth_eigenvector, label='Fifth Eigenvector', marker='o')
plt.plot(x_range, sin_function, label=r'$\sqrt{2} \sin(\pi x)$', linestyle='--')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Fifth Eigenvector and Sin Function for 10x10 Matrix')
plt.grid(True)
plt.show()
