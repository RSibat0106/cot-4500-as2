# Problem 1: Neville's Method for f(3.7)
def neville_method(x_values, y_values, x_interp):
    n = len(x_values)
    Q = [[0.0] * n for _ in range(n)]

    for i in range(n):
        Q[i][0] = y_values[i]

    for j in range(1, n):
        for i in range(n - j):
            Q[i][j] = ((x_interp - x_values[i + j]) * Q[i][j - 1] +
                       (x_values[i] - x_interp) * Q[i + 1][j - 1]) / (x_values[i] - x_values[i + j])

    return Q[0][n - 1]

x_values_neville = [3.6, 3.8, 3.9]
y_values_neville = [1.675, 1.436, 1.318]
x_interp_neville = 3.7

result_neville = neville_method(x_values_neville, y_values_neville, x_interp_neville)
print(result_neville)
print("")
# Problem 2: Newton's Forward Difference Table
def newton_forward_difference_table(x_values, y_values):
    n = len(y_values)
    diff_table = [[0] * n for _ in range(n)]
    for i in range(n):
        diff_table[i][0] = y_values[i]

    for j in range(1, n):
        for i in range(n - j):
            diff_table[i][j] = diff_table[i + 1][j - 1] - diff_table[i][j - 1]

    return diff_table

def newton_forward_polynomial(x_values, diff_table):
    h = x_values[1] - x_values[0]
    f0 = diff_table[0][0]
    f1 = diff_table[0][1] / h
    f2 = diff_table[0][2] / (2 * h ** 2)
    f3 = diff_table[0][3] / (6 * h ** 3)
    return f0, f1, f2, f3

x_values_newton = [7.2, 7.4, 7.5, 7.6]
y_values_newton = [23.5492, 25.3913, 26.8224, 27.4589]

difference_table_newton = newton_forward_difference_table(x_values_newton, y_values_newton)
coefficients_newton = newton_forward_polynomial(x_values_newton, difference_table_newton)

for coef in coefficients_newton:
    print(coef)
    print("")

# Problem 3: Approximate f(7.3) using Newton’s Forward Method
def newton_forward_interpolation(x_values, diff_table, x_interp):
    h = x_values[1] - x_values[0]
    p = (x_interp - x_values[0]) / h
    f_x = diff_table[0][0]
    p_term = 1

    for j in range(1, len(x_values)):
        p_term *= (p - (j - 1)) / j
        f_x += p_term * diff_table[0][j]

    return f_x

x_interp_newton = 7.3
result_newton_interp = newton_forward_interpolation(x_values_newton, difference_table_newton, x_interp_newton)
print(result_newton_interp)
print("")

# Problem 4: Hermite Divided Difference Table
def hermite_divided_difference(x_values, y_values, dy_values):
    n = len(x_values)
    size = 2 * n
    H = [[0] * (size + 1) for _ in range(size)]
    z = [x for x in x_values for _ in (0, 1)]

    for i in range(n):
        H[2 * i][0] = H[2 * i + 1][0] = y_values[i]
        H[2 * i + 1][1] = dy_values[i]
        if i != n - 1:
            H[2 * i][1] = (y_values[i + 1] - y_values[i]) / (x_values[i + 1] - x_values[i])

    for j in range(2, size):
        for i in range(size - j):
            H[i][j] = (H[i + 1][j - 1] - H[i][j - 1]) / (z[i + j] - z[i])

    return H

x_values_hermite = [3.6, 3.8, 3.9]
y_values_hermite = [1.675, 1.436, 1.318]
dy_values_hermite = [-1.195, -1.188, -1.182]

hermite_table = hermite_divided_difference(x_values_hermite, y_values_hermite, dy_values_hermite)
for row in hermite_table:
    print(row)
    print("")

# Problem 5: Cubic Spline Interpolation - Find Matrices A, b, x
def cubic_spline_matrices(x_values, y_values):
    n = len(x_values) - 1
    A = [[0] * (n + 1) for _ in range(n + 1)]
    b = [0] * (n + 1)

    for i in range(1, n):
        A[i][i - 1] = 1
        A[i][i] = 4
        A[i][i + 1] = 1
        b[i] = 3 * ((y_values[i + 1] - y_values[i]) / (x_values[i + 1] - x_values[i]) -
                    (y_values[i] - y_values[i - 1]) / (x_values[i] - x_values[i - 1]))

    A[0][0] = 1
    A[n][n] = 1

    # Solve Ax = b using Gaussian elimination
    x_vector = [0] * (n + 1)
    for i in range(n + 1):
        for j in range(i + 1, n + 1):
            if A[i][i] != 0:
                factor = A[j][i] / A[i][i]
                for k in range(n + 1):
                    A[j][k] -= factor * A[i][k]
                b[j] -= factor * b[i]

    for i in range(n, -1, -1):
        sum_ax = sum(A[i][j] * x_vector[j] for j in range(i + 1, n + 1))
        if A[i][i] != 0:
            x_vector[i] = (b[i] - sum_ax) / A[i][i]

    return A, b, x_vector

x_values_spline = [2, 5, 8, 10]
y_values_spline = [3, 5, 7, 9]

matrix_A, vector_b, vector_x = cubic_spline_matrices(x_values_spline, y_values_spline)

#Matrix M
for row in matrix_A:
    print(row)
    print("")

#Vector B
print(vector_b)
print("")
#Vector X
print(vector_x)
print("")
