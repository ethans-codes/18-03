import numpy as np
from colors import bcolors
import sympy as sp

def max_steps(a, b, err):
    """
    Receives 3 parameters:
        1.  a - start value.
        2.  b - end  value.
        3.  err - value of tolerable error

    Returns variables:
        1.  S - The minimum number of iterations required to reach the desired accuracy
    """
    s = int(np.ceil(- np.log2(err / (b - a)) / np.log2(2) - 1))
    return s

def find_intervals_between(f, a, b, step = 0.1):
    intervals = []
    cur = a
    while (cur <= b):
        if f(cur) * f(cur + step) <0:
            intervals.append((cur, cur + step))
        cur += step

    return intervals

def Der(f):
    x = sp.symbols('x')
    # Convert the function to a SymPy expression
    f_expr = sp.sympify(f(x))

    # Compute the derivative
    f_prime = sp.diff(f_expr, x)

    # Convert f_prime to a Python function
    f_prime_func = sp.lambdify(x, f_prime, 'numpy')

    return f_prime_func

def Run_Bisection_Method(f, a, b, tol=1e-6):

    if np.sign(f(a)) != np.sign(f(b)):
        roots = bisection_method(f, a, b, tol)
        print(bcolors.OKBLUE, f"\nThe equation f(x) has an approximate root at x = {roots}", bcolors.ENDC)

    else:
        pos_roots = find_intervals_between(f, a, b)
        for root in pos_roots:
            a, b = root
            tmp = bisection_method(f, a, b, tol)
            if tmp != None:
                print(bcolors.OKBLUE, f"\nThe equation f(x) has an approximate root at x = {tmp}", bcolors.ENDC)
            else:
                print(bcolors.OKBLUE, "No root found in the interval.", bcolors.ENDC)

        # Calculation for 𝑓′(𝑥) - NOT FINISHED!!!
        pos_ft_roots = []
        for current_range, next_range in zip(pos_roots, pos_roots[1:]):
            tmp = (current_range[1], next_range[0])  # Create the new range
            pos_ft_roots.append(tmp)
        #f_t=Der(f)


def bisection_method(f, a, b, tol=1e-6):
    """
    Performs Iterative methods for Nonlinear Systems of Equations to determine the roots of the given function f
    Receives 4 parameters:
        1. f - continuous function on the interval [a, b], where f (a) and f (b) have opposite signs.
        2. a - start value.
        3. b - end  value.
        4. tol - the tolerable error , the default value will set as 1e-16

    Returns variables:
        1.  c - The approximate root of the function f
    """

    c, k = 0, 0
    steps = max_steps(a, b, tol)  # calculate the max steps possible
    print("{:<10} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}".format("Iteration", "a", "b", "f(a)", "f(b)", "c", "f(c)"))

    while abs(b - a) >= tol and k < steps:  # Adjusted stopping condition
        c = a + (b - a) / 2
        if f(c) == 0:
            return c  # Procedure completed successfully

        # if sign changed between steps
        if (f(c) * f(a)) < 0:
            # move forward
            b = c
        else:
            # move backward
            a = c

        print("{:<10} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f} {:<15.6f}".format(k, a, b, f(a), f(b), c, f(c)))
        if round(a) == b:
            return None
        k += 1
    return c  # return the current root

if __name__ == '__main__':
    #f = lambda x: x ** 4 + x**3 - 3*x**2
    f = lambda x: 6 * x ** 4 - 7 * x ** 3 - 2 * x + 1
    Run_Bisection_Method(f, 0, 5, 10**-10)