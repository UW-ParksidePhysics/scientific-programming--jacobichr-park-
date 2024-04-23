import numpy as np
import matplotlib.pyplot as plt


# Define the matrix H using NumPy diagflat and ones functions and matrix operations
s = 1 / (2 * (1 / (5 + 1))**2)
values = np.array([2 - 1000, -12 - 100, 0 - 12 - 10, 0 - 12 - 1, 0 - 12])
H = s * np.diagflat(values)

# Find the eigenvalues and eigenvectors of H using NumPy linalg eig function
eigenvalues, eigenvectors = np.linalg.eig(H)

# Plot the fifth eigenvector as the y-values against
# x_values = np.linspace(1/(5+1), 5/(5+1), 5); together with the function
# √2 sin(πx) on the range 0 < x <= 1.
x_values = np.linspace(1 / (5 + 1), 5 / (5 + 1), 5)
fifth_eigenvector = eigenvectors[:, 4] 

# Define the range for the sin function plot
x_range = np.linspace(0, 1, 100)
sin_function = np.sqrt(2) * np.sin(np.pi * x_range)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_values, fifth_eigenvector, label='Fifth Eigenvector', marker='o')
plt.plot(x_range, sin_function, label=r'$\sqrt{2} \sin(\pi x)$', linestyle='--')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Fifth Eigenvector and Sin Function')
plt.grid(True)
plt.show()
