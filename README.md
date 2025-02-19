# cot-4500-as2
Programming Assignment 2

# Numerical Methods Programming Assignment

## **Overview**
This repository contains an implementation of various **numerical interpolation techniques** used to approximate function values based on given data points. The following methods are implemented:

- **Neville’s Method** for polynomial interpolation
- **Newton’s Forward Method** for polynomial approximations
- **Hermite Polynomial Approximation** using the Divided Difference Method
- **Cubic Spline Interpolation** to construct smooth curves

This assignment is designed to work with **Python IDLE 3.13.3**, ensuring compatibility without external libraries like NumPy.

---

## **Methods Implemented**
### **1. Neville’s Method**
- Uses an iterative approach to approximate `f(3.7)`.
- Constructs **Neville’s interpolation table** to compute results.

### **2. Newton’s Forward Method**
- Constructs the **divided difference table**.
- Computes polynomial approximations for **degrees 1, 2, and 3**.
- Approximates `f(7.3)` using **Newton’s Forward Interpolation**.

### **3. Hermite Polynomial Approximation**
- Implements the **divided difference method** for Hermite interpolation.
- Computes the **Hermite divided difference table** for given `(x, f(x), f'(x))` values.

### **4. Cubic Spline Interpolation**
- Constructs **Matrix A** and **Vector b**.
- Solves for **Vector x** using **Gaussian Elimination**.
- Ensures continuity and smoothness between data points.

---

## **Installation**
This program is implemented in **pure Python (without NumPy)**, making it compatible with **Python 3.13.3**. 

To install Python 3.13.3 (if not already installed):

1. Download Python 3.13.3 from the [official Python website](https://www.python.org/downloads/).
2. Follow the installation instructions for your operating system.

---

## **How to Run the Program**
You can run the script **directly from the command line** using the following steps:

### **1. Clone the Repository**
```sh
git clone https://github.com/yourusername/numerical-methods-assignment.git
cd numerical-methods-assignment

run assignment_2.py
