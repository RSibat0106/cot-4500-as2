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

def newton_forward_difference_table(x_values, y_values):
    """
    Constructs the forward difference table for Newton's forward interpolation.
    
    :param x_values: List of x values (equally spaced)
    :param y_values: List of corresponding f(x) values
    :return: Forward difference table
    """
    n = len(y_values)
    diff_table = [y_values[:]]  # First column is f(x) values

    for j in range(1, n):
        column = []
        for i in range(n - j):
            diff = diff_table[j - 1][i + 1] - diff_table[j - 1][i]
            column.append(diff)
        diff_table.append(column)
    
    return diff_table

def print_forward_difference_table(diff_table):
    """
    Prints the forward difference table in a readable format.
    """
    print("\nNewton's Forward Difference Table:")
    for row in diff_table:
        print([round(val, 6) for val in row])

def newton_forward_polynomial(x_values, diff_table):
    """
    Constructs and prints Newton's forward interpolation polynomials of degree 1, 2, and 3.
    
    :param x_values: List of x values (equally spaced)
    :param diff_table: Forward difference table
    """
    h = x_values[1] - x_values[0]  # Step size
    x = x_values[0]  # Starting x value
    
    # Compute coefficients
    f0 = diff_table[0][0]
    f1 = diff_table[1][0] / h
    f2 = diff_table[2][0] / (2 * h**2)
    f3 = diff_table[3][0] / (6 * h**3)
    
    # Print polynomial approximations
    print("\nNewton's Forward Interpolation Polynomials:")
    
    # Degree 1 Polynomial: f(x) ≈ f0 + f1*(x - x0)
    print(f"P1(x) = {f0:.6f} + ({f1:.6f}) * (x - {x:.1f})")
    
    # Degree 2 Polynomial: P1(x) + f2*(x - x0)*(x - x1)
    print(f"P2(x) = {f0:.6f} + ({f1:.6f}) * (x - {x:.1f}) + ({f2:.6f}) * (x - {x:.1f}) * (x - {x_values[1]:.1f})")
    
    # Degree 3 Polynomial: P2(x) + f3*(x - x0)*(x - x1)*(x - x2)
    print(f"P3(x) = {f0:.6f} + ({f1:.6f}) * (x - {x:.1f}) + ({f2:.6f}) * (x - {x:.1f}) * (x - {x_values[1]:.1f})"
          f" + ({f3:.6f}) * (x - {x:.1f}) * (x - {x_values[1]:.1f}) * (x - {x_values[2]:.1f})")

# Given data points
x_values = [7.2, 7.4, 7.5, 7.6]
y_values = [23.5492, 25.3913, 26.8224, 27.4589]

# Compute forward difference table
difference_table = newton_forward_difference_table(x_values, y_values)

# Print the table
print_forward_difference_table(difference_table)

# Compute and print polynomials of degree 1, 2, and 3
newton_forward_polynomial(x_values, difference_table)
