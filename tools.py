import numpy as np
import matplotlib.pyplot as plt

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

# Find maximum


def maximum(tab):
    idx_max = 0
    max = tab[0]
    for i in range(0, len(tab)):
        if tab[i] > max:
            max = tab[i]
            idx_max = i
    return max, idx_max

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

# Schottky


def schottky(T, E, n=1, k=1.380649e-23):
    """Calculate the Schottky anomaly"""
    """T: temperature in Kelvin, E: energy in Joules, n: number of particles, k: Boltzmann constant"""
    x = (E)/(k*T)
    cs = k*(x**2)*(np.exp(x)/(1+np.exp(x))**2)
    return n*cs


def dev_schottky(T, E, n=1, k=1.380649e-23):
    """Calculate the derivative of the Schottky anomaly"""
    a = E / k
    exp_at = np.exp(a / T)
    num = (2 * T + a) * exp_at + (2 * T - a) * (exp_at ** 2)
    denom = (T ** 4) * (1 + exp_at) ** 3
    return - n * k * (a ** 2) * num / denom


def max_schottky(x, y, min, max):
    """Find the maximum of the Schottky anomaly in a given interval"""
    x_interval, y_interval = tab_interval(x, y, min, max)
    maxi, i = maximum(y_interval)
    x_maxi = x_interval[i]
    return x_maxi, maxi

