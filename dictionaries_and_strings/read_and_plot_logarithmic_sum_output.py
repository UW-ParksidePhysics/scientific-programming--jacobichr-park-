import re
import matplotlib.pyplot as plt

# Function to parse the output file
def parse_sum_output(filepath):
    tolerances, errors, maximum_indices = [], [], []
    with open(filepath, 'r') as file:
        for line in file:
            # Use regular expression to extract numbers
            match = re.search(r'epsilon: (\d+e-\d+), exact error: (\d+.\d+e-\d+), n=(\d+)', line)
            if match:
                epsilon, error, n = match.groups()
                tolerances.append(float(epsilon))
                errors.append(float(error))
                maximum_indices.append(int(n))
    return tolerances, errors, maximum_indices

# Function to plot the error
def plot_logarithmic_sum_error(tolerances, errors, maximum_indices):
    plt.figure()
    plt.semilogy(maximum_indices, errors, marker='o')
    plt.xlabel('Maximum Index (n)')
    plt.ylabel('Approximation Error (Δ)')
    plt.title('Logarithmic Plot of Approximation Error vs Maximum Index')
    plt.grid(True)
    plt.show()

# Test the function
if __name__ == '__main__':
    tolerances, errors, maximum_indices = parse_sum_output('dictionaries_and_strings/logarithmic_sum.out')
    plot_logarithmic_sum_error(tolerances, errors, maximum_indices)
