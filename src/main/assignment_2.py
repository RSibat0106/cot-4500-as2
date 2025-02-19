def neville_method(x_data, y_data, x_interp):
    """
    Implements Neville's interpolation method to find f(x_interp).
    
    :param x_data: List of x values (known data points)
    :param y_data: List of corresponding f(x) values
    :param x_interp: The x value to interpolate
    :return: Interpolated value at x_interp
    """
    n = len(x_data)
    Q = [[0.0] * n for _ in range(n)]  # Initialize empty table

    # Fill the first column with function values
    for i in range(n):
        Q[i][0] = y_data[i]

    # Apply Neville's formula
    for j in range(1, n):
        for i in range(n - j):
            Q[i][j] = ((x_interp - x_data[i + j]) * Q[i][j - 1] + 
                       (x_data[i] - x_interp) * Q[i + 1][j - 1]) / (x_data[i] - x_data[i + j])

    return Q[0][n - 1]  # Top-right value is the interpolated result

# Given data points
x_values = [3.6, 3.8, 3.9]
y_values = [1.675, 1.436, 1.318]

# Interpolation point
x_interp = 3.7

# Compute interpolated value
result = neville_method(x_values, y_values, x_interp)
print(f"Interpolated value at x = {x_interp}: {result:.6f}")

import numpy as np
import pandas as pd

def newton_forward_difference_table(x_values, y_values):
    """
    Constructs the forward difference table for Newton's forward interpolation (pure Python).
    
    :param x_values: List of x values (equally spaced)
    :param y_values: List of corresponding f(x) values
    :return: Forward difference table as a list of lists
    """
    n = len(y_values)
    diff_table = [[0] * n for _ in range(n)]  # Create an empty table
    for i in range(n):
        diff_table[i][0] = y_values[i]  # Fill first column with f(x) values

    # Compute forward differences
    for j in range(1, n):
        for i in range(n - j):
            diff_table[i][j] = diff_table[i + 1][j - 1] - diff_table[i][j - 1]

    return diff_table

def print_difference_table(x_values, diff_table):
    """
    Prints the forward difference table in a readable format.
    """
    print("\nNewton's Forward Difference Table:")
    n = len(x_values)
    headers = ["x", "f(x)"] + [f"Δ^{j} f(x)" for j in range(1, n)]
    print("{:<8} {:<12}".format(headers[0], headers[1]), end="")
    for h in headers[2:]:
        print("{:<12}".format(h), end="")
    print()

    for i in range(n):
        print(f"{x_values[i]:<8} {diff_table[i][0]:<12.6f}", end="")
        for j in range(1, n - i):
            print(f"{diff_table[i][j]:<12.6f}", end="")
        print()

def newton_forward_polynomial(x_values, diff_table):
    """
    Constructs and prints Newton's forward interpolation polynomials of degree 1, 2, and 3.
    
    :param x_values: List of x values (equally spaced)
    :param diff_table: Forward difference table
    """
    h = x_values[1] - x_values[0]  # Step size
    x0 = x_values[0]  # First x value
    
    # Extract first-row forward differences
    f0 = diff_table[0][0]
    f1 = diff_table[0][1] / h
    f2 = diff_table[0][2] / (2 * h**2)
    f3 = diff_table[0][3] / (6 * h**3)
    
    # Create polynomial strings
    P1 = f"P1(x) = {f0:.6f} + ({f1:.6f}) * (x - {x0:.1f})"
    P2 = f"P2(x) = {f0:.6f} + ({f1:.6f}) * (x - {x0:.1f}) + ({f2:.6f}) * (x - {x0:.1f}) * (x - {x_values[1]:.1f})"
    P3 = (f"P3(x) = {f0:.6f} + ({f1:.6f}) * (x - {x0:.1f}) + ({f2:.6f}) * (x - {x0:.1f}) * (x - {x_values[1]:.1f})"
          f" + ({f3:.6f}) * (x - {x0:.1f}) * (x - {x_values[1]:.1f}) * (x - {x_values[2]:.1f})")
    
    return [P1, P2, P3]

# Given data points
x_values = [7.2, 7.4, 7.5, 7.6]
y_values = [23.5492, 25.3913, 26.8224, 27.4589]

# Compute forward difference table
difference_table = newton_forward_difference_table(x_values, y_values)

# Print the difference table
print_difference_table(x_values, difference_table)

# Compute polynomial approximations
polynomials = newton_forward_polynomial(x_values, difference_table)

# Display polynomial approximations
print("\nNewton's Forward Interpolation Polynomials:")
for poly in polynomials:
    print(poly)
