import sympy as sp

def construct_polynomial_from_roots(x, roots):
    # Initialize the polynomial to 1 (neutral element for multiplication)
    p = 1
    # Loop through all roots and construct the polynomial
    for root in roots:
        p *= (x - root)
    return p

if __name__ == "__main__":
  # Example of how to use the function
  x = 42  # Define the symbol for x
  roots = [1, 2, 3, 7, 27939]  # Example roots of the polynomial
  polynomial = construct_polynomial_from_roots(x, roots)

  print(polynomial)