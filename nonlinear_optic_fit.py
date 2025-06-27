import numpy as np
import matplotlib.pyplot as plt
import tools
import constants as cnt
import scipy.optimize as opt
from schottky_analysis import schottky

# Fit function


def model_optic(x, beta, gamma, n, E, nu):
    """Fit function for the nonlinear analysis of C/T - C_schottky vs T².
    Parameters: x = T² (K²), beta = mJ/K⁴.mol, gamma = mJ/K².mol, n = dimensionless"""
    phonon = beta * x + nu*(x**2)
    y = np.sqrt(x)
    # schottky = (y**2) * np.exp(y)/((np.exp(y) + 1)**2)
    cs = schottky(y, E, n)/y
    return phonon + gamma + n*cnt.r*1e3*cs

# Non linear fit


def optic_fit(a, b, x_carre, y, err_y, bounds):
    """Non linear fit usinf curve fit from scipy library : 
    a, b = sclars (bounds), 
    x_carre = array
    y = array
    err_y = array
    bounds = bounds (2-tuple of arrays-like)"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    fit = opt.curve_fit(model_optic, x_carre_interval, y_interval, bounds=bounds,
                        absolute_sigma=True)
    return fit[0]


def plot_fit_optic(a, b, x_carre, y, err_y, bounds):
    """Plotting data and non linear fit usinf curve fit from scipy library : 
    a, b = sclars (bounds in squared temperature), 
    x_carre = array
    y = array
    err_y = array
    bounds = bounds (2-tuple of arrays-like)"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    beta, gamma, n, E, nu = optic_fit(
        a, b, x_carre, y, err_y, bounds=bounds)
    print("Beta, Gamma, n, E, nu : ", beta, gamma, n, E, nu)
    plt.figure()
    plt.plot(x_carre_interval, y_interval, "g.", label="Experimental")
    plt.plot(x_carre_interval, model_optic(
        x_carre_interval, beta, gamma, n, E, nu), "c-", label="Fit")
    plt.xlabel("Temperature (K)")
    plt.ylabel("C/T (mJ/K²/mol)")
    plt.grid(True)
    plt.legend()
    plt.show()


def main():
    pass


if __name__ == "__main__":
    main()
