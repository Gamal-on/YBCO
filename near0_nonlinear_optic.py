import numpy as np
import matplotlib.pyplot as plt
import tools
import constants as cnt
import scipy.optimize as opt

# Fit function


def model_optic(x, beta, gamma, n, E, nu):
    """Fit function for the nonlinear analysis of C/T - C_schottky vs T².
    Parameters: x = T² (K²), beta = mJ/K⁴.mol, gamma = mJ/K².mol, n = dimensionless"""
    phonon = beta * x + nu*(x**2)
    y = E/(cnt.k*np.sqrt(x))
    schottky = (y**2) * np.exp(y)/((np.exp(y) + 1)**2)
    return phonon + gamma + n*cnt.r*1e3*schottky/np.sqrt(x)

# Non linear fit


def optic_fit(a, b, x_carre, y, err_y, bounds=([0.1, 0, 5e-3, 9.3e-23], [1, 5, 1.1e-2, 1.2e-22])):
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


def plot_fit_optic(a, b, x_carre, y, err_y, bounds=([0.1, 0, 5e-3, 9.3e-23], [1, 5, 1.1e-2, 1.2e-22])):
    """Plotting data and non linear fit usinf curve fit from scipy library : 
    a, b = sclars (bounds), 
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
    plt.grid(True)
    plt.legend()
    plt.show()


def main():
    plot_fit_optic(0, 400, cnt.squared_temperature, cnt.C_div_T, cnt.err_C_divT,
                   bounds=([0.1, 0, 5e-3, 9.8e-23, 0], [1, 5, 1.1e-2, 1.2e-22, 1]))


if __name__ == "__main__":
    main()
