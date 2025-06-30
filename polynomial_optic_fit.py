import numpy as np
import matplotlib.pyplot as plt
import tools
import constants as cnt
import scipy.optimize as opt
from schottky_analysis import schottky

# Fit function


def model_quadra(x, alpha, beta, gamma):
    """Fit function for the nonlinear analysis of C/T - C_schottky vs T².
    Parameters: x = T² (K²), beta = mJ/K⁴.mol, gamma = mJ/K².mol, n = dimensionless"""
    phonon = beta * x
    quadra = alpha*(x**2)
    return phonon + gamma + quadra

# Non linear fit


def nonlinear_fit_quadra_opt(a, b, x_carre, y, bounds):
    """Non linear fit using curve fit from scipy library :
    a, b = sclars (bounds), WARNING : squared bounds
    x_carre = array
    y = array
    err_y = array
    bounds = bounds (2-tuple of arrays-like), alpha beta gamma"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    fit = opt.curve_fit(model_quadra, x_carre_interval, y_interval, bounds=bounds,
                        absolute_sigma=True)
    return fit[0]


def nonlinear_fit_quadra_pol(a, b, x_carre, y):
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    fit = np.polyfit(x_carre_interval, y_interval, deg=2)
    return fit


def plot_fit_quadra(a, b, x_carre, y, bounds, opt=True):
    """Plotting data and non linear fit usinf curve fit from scipy library :
    a, b = sclars (bounds), WARNING : squared bounds
    x_carre = array
    y = array
    err_y = array
    bounds = bounds (2-tuple of arrays-like)
    opt : bool"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    if opt:
        alpha, beta, gamma = nonlinear_fit_quadra_opt(
            a, b, x_carre_interval, y_interval, bounds)
    else:
        alpha, beta, gamma = nonlinear_fit_quadra_pol(
            a, b, x_carre_interval, y_interval)
    print("Alpha, Beta, Gamma : ", alpha, beta, gamma)
    plt.figure()
    plt.plot(x_carre_interval, y_interval, "g.", label="Experimental")
    plt.plot(x_carre_interval, model_quadra(
        x_carre_interval, alpha, beta, gamma), "c-", label="Fit")
    plt.grid(True)
    plt.legend()
    plt.show()


def main():
    plot_fit_quadra(144, 400, cnt.squared_temperature_HPHT,
                    cnt.hc_div_temp_HPHT, bounds=([0, 0.1, 0], [1, 1, 10]), opt=False)


if __name__ == "__main__":
    main()
