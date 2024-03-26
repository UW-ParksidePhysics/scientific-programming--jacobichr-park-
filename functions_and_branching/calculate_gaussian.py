import math

def gaussian(x, mu, sigma):
    # Ensure sigma is positive to avoid division by zero or negative standard deviation
    if sigma <= 0:
        raise ValueError("Standard deviation must be positive")

    # Compute the Gaussian function
    factor = 1 / (sigma * math.sqrt(2 * math.pi))
    exponent = math.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
    return factor * exponent

def print_gaussian_table(mu, sigma):
  start = mu - sigma
  end = mu + sigma
  step = sigma / 10  # Define a reasonable step based on the standard deviation

  # Print the header of the table
  print(f"mu={mu} sigma={sigma}")
  print(f"{'x':>10} | {'Gaussian(x)':>20}")
  print("-" * 33)

  # Calculate and print the Gaussian function values for the range
  x = start
  while x <= end:
      result = gaussian(x, mu, sigma)
      print(f"{x:>10.2f} | {result:>20.6f}")
      x += step

if __name__ == "__main__":
  # Example usage
  mu = 42  # Mean
  sigma = 21  # Standard deviation
  print_gaussian_table(mu, sigma)