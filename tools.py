import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import constants as cnt


# Function to choose a certain interval

def tab_interval(x, y, min, max):
    """Choose bounds in x, return the corresponding values in x and y"""
    idx_min = 0
    idx_max = 0
    tab = []
    for valeur in x:
        if valeur < min:
            idx_min += 1
            idx_max += 1
        elif valeur == min:
            tab.append(valeur)
            idx_max += 1
        elif min < valeur < max:
            tab.append(valeur)
            idx_max += 1
        elif valeur == max:
            tab.append(valeur)
            idx_max += 1
    y_interval = y[idx_min: idx_max]
    x_interval = np.array(tab)
    return x_interval, y_interval

# Find extremum


def maximum(tab):
    idx_max = 0
    max = tab[0]
    for i in range(0, len(tab)):
        if tab[i] > max:
            max = tab[i]
            idx_max = i
    return max, idx_max


def minimum(tab):
    idx_min = 0
    min = tab[0]
    for i in range(0, len(tab)):
        if tab[i] < min:
            min = tab[i]
            idx_min = i
    return min, idx_min

# Plot


def plot(x, y, min, max):
    x_interval, y_interval = tab_interval(x, y, min, max)
    plt.figure()
    plt.plot(x_interval, y_interval, ".g")
    plt.grid(True)
    plt.show()

# Newton Raphson


def resol_nr(f, fp, N, xstart, eps):
    """
    Solve f(x) = 0 using Newton-Raphson method
    f: function to solve
    fp: derivative of f
    N: maximum number of iterations
    xstart: initial guess
    eps: tolerance for convergence"""
    i = 0
    x2 = xstart
    condition = True
    while condition:
        i += 1
        x1 = x2
        x2 = x2-f(x2)/fp(x2)
        condition = np.abs(x2-x1) > eps and N < i
    return x2, f(x2)

# Dichotomie


def dicho_1(f, eps, a, b):
    """
    Solve f(x) = 0 using the dichotomy method
    f: function to solve
    eps: tolerance for convergence
    a: lower bound of the interval
    b: upper bound of the interval"""
    borne_min = a
    borne_max = b
    milieu = (b-a)/2
    res = f(milieu)
    while np.abs(res) > eps:
        resmin = f(borne_min)
        if resmin * res > 0:
            borne_min = milieu
        else:
            borne_max = milieu
        milieu = (borne_min + borne_max)/2
        res = f(milieu)
    return res, milieu

# Integrate


def simpson(f, a, b, n):
    somme = 0
    dx = (b-a)/n
    for i in range(0, n):
        x = a + i*dx
        x1 = x + dx
        somme += (f(x) + 4*f((x + x1)/2) + f(x1))
    return (b-a)/(6*n) * somme


def integrate(x, y, a, b):
    """Integrate y with respect to x from a to b using the trapezoidal rule."""
    mask = (x >= a) & (x <= b)
    x_masked = x[mask]
    y_masked = y[mask]
    integral = cumulative_trapezoid(y_masked, x_masked, initial=0)
    return integral

# Debye temperature


def debye_temperature(beta):
    """
    Calculate the Debye temperature and gamma from the linear fit parameters
    Returns the Debye temperature in K, gamma in J/K².mol and their respectiv errors"""
    pi4 = np.pi**4
    factor = cnt.N*12*cnt.k*pi4/5
    temp_debye = np.cbrt(factor/(beta*1e-3))
    return temp_debye

def main():
    print(debye_temperature(0.7664))


if __name__ == "__main__":
    main()

