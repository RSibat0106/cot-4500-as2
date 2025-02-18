python3 -m pip install numpy
import numpy as np

def neville_interpolation(x_data, y_data, x_interp):
    """
    Neville's method for polynomial interpolation.
    
    :param x_data: List of x values
    :param y_data: List of corresponding f(x) values
    :param x_interp: The x value to interpolate
    :return: Interpolated value at x_interp
    """
    n = len(x_data)
    Q = np.zeros((n, n))
    
    # Initialize first column with function values
    for i in range(n):
        Q[i, 0] = y_data[i]

    # Compute Neville’s table
    for j in range(1, n):
        for i in range(n - j):
            Q[i, j] = ((x_interp - x_data[i + j]) * Q[i, j - 1] + 
                       (x_data[i] - x_interp) * Q[i + 1, j - 1]) / (x_data[i] - x_data[i + j])
    
    return Q[0, n - 1]  # Return the top-right value of the table

# Given data points
x_values = [3.6, 3.8, 3.9]
y_values = [1.675, 1.436, 1.318]

# Interpolation point
x_interp = 3.7

# Compute interpolation
interpolated_value = neville_interpolation(x_values, y_values, x_interp)
print(f"Interpolated value at x = {x_interp}: {interpolated_value:.6f}")
