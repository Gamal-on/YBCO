import numpy as np
import matplotlib.pyplot as plt
import tools
import constants as cnt
import scipy.optimize as opt
from schottky_analysis import schottky

# Fit function


def model_quadra(x, beta, gamma, alpha):
    """Fit function for the nonlinear analysis of C/T - C_schottky vs T².
    Parameters: x = T² (K²), beta = mJ/K⁴.mol, gamma = mJ/K².mol, n = dimensionless"""
    phonon = beta * x
    quadra = alpha*(x**2)
    return phonon + gamma + quadra

# Non linear fit


def nonlinear_fit_quadra(a, b, x_carre, y, bounds):
    """Non linear fit usinf curve fit from scipy library :
    a, b = sclars (bounds),
    x_carre = array
    y = array
    err_y = array
    bounds = bounds (2-tuple of arrays-like)"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    fit = opt.curve_fit(model_quadra, x_carre_interval, y_interval, bounds=bounds,
                        absolute_sigma=True)
    return fit


def plot_fit_quadra(a, b, x_carre, y, bounds):
    """Plotting data and non linear fit usinf curve fit from scipy library :
    a, b = sclars (bounds),
    x_carre = array
    y = array
    err_y = array
    bounds = bounds (2-tuple of arrays-like)"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    beta, gamma, alpha = nonlinear_fit_quadra(
        a, b, x_carre, y, bounds)[0]
    print("Beta, Gamma, alpha : ", beta, gamma, alpha)
    plt.figure()
    plt.plot(x_carre_interval, y_interval, "g.", label="Experimental")
    plt.plot(x_carre_interval, model_quadra(
        x_carre_interval, beta, gamma, alpha), "c-", label="Fit")
    plt.grid(True)
    plt.legend()
    plt.show()


def main():
    pass


if __name__ == "__main__":
    main()
