import numpy as np
import matplotlib.pyplot as plt
import tools

# Constantes et tableaux

k = 1.380649e-23
delta = 2.9461005*k
temp = np.arange(0, 3, 1e-3)


def max_schottky(x, y, min, max):
    """Find the maximum of the Schottky anomaly in a given interval"""
    x_interval, y_interval = tools.tab_interval(x, y, min, max)
    maxi, i = tools.maximum(y_interval)
    x_maxi = x_interval[i]
    return x_maxi, maxi


def schottky(T, E, n=1, r=8.31446261815324, na=6.02214076e23, k=1.380649e-23):
    """Calculate the Schottky anomaly"""
    """T: temperature in Kelvin, E: energy in Joules, n: number of particles, k: Boltzmann constant"""
    x = (E)/(k*T)
    cs = (x**2)*(np.exp(x)/(1+np.exp(x))**2)
    return r*cs


def dev_schottky(T, E, n=1):
    a = E / k
    exp_at = np.exp(a / T)
    num = (2 * T + a) * exp_at + (2 * T - a) * (exp_at ** 2)
    denom = (T ** 4) * (1 + exp_at) ** 3
    return - n * k * (a ** 2) * num / denom


def plot_schottky(T, E, n=1):
    plt.figure()
    plt.plot(T, schottky(T, E, n=1), "r")
    plt.grid(True)
    plt.show()


def plot_dev_schottky(T, E, n=1):
    plt.figure()
    plt.plot(T, dev_schottky(T, E, n=1), "r")
    plt.grid(True)
    plt.show()

# Relation entre T* et E en recherchant les zeros de la dérivée : E = k*T*2.399357318878174, méthode dichotomie


def f(x):
    return 2+x+(2-x)*np.exp(x)


def fp(x):
    return 1 + 2*np.exp(x) - np.exp(x) - x*np.exp(x)
