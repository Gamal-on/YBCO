import numpy as np
import matplotlib.pyplot as plt
import tools
import schottky_analysis
import constants as cnt
import scipy.optimize as opt

# Fit function


def fit_func(x, beta, gamma, n, E):
    """Fit function for the nonlinear analysis of C/T - C_schottky vs T².
    Parameters: x = T² (K²), beta = mJ/K⁴.mol, gamma = mJ/K².mol, n = dimensionless"""
    phonon = beta * x
    y = E/(cnt.k*np.sqrt(x))
    schottky = (y**2) * np.exp(y)/((np.exp(y) + 1)**2)
    return phonon + gamma + n*cnt.r*1e3*schottky/np.sqrt(x)

# Non linear fit


def nonlinear_fit(a, b, x_carre, y, err_y, bounds=([0.1, 0, 1e-3, 9e-23], [1, 10, 1e-2, 1.2e-22])):
    """Non linear fit usinf curve fit from scipy library : 
    a, b = sclars (bounds), 
    x_carre = array
    y = array
    err_y = array
    bounds = bounds (2-tuple of arrays-like)"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    fit = opt.curve_fit(fit_func, x_carre_interval, y_interval, bounds=bounds,
                        sigma=err_y, absolute_sigma=True)
    return fit[0]


def plot_fit(a, b, x_carre, y, err_y, bounds=([0.1, 0, 1e-3, 9e-23], [1, 10, 1e-2, 1.2e-22])):
    """Plotting data and non linear fit usinf curve fit from scipy library : 
    a, b = sclars (bounds), 
    x_carre = array
    y = array
    err_y = array
    bounds = bounds (2-tuple of arrays-like)"""
    x_carre_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    beta, gamma, n, E = nonlinear_fit(a, b, x_carre, y, err_y, bounds=bounds)
    print("Beta, Gamma, n, E : ", beta, gamma, n, E)
    plt.figure()
    plt.plot(x_carre_interval, y_interval, "g.", label="Experimental")
    plt.plot(x_carre_interval, fit_func(
        x_carre_interval, beta, gamma, n, E), "c-", label="Fit")
    plt.grid(True)
    plt.legend()
    plt.show()

# Debye temperature


def debye_temperature(a, b, N=78e23):
    """
    Calculate the Debye temperature from the fit parameters
    Returns the Debye temperature in K and gamma in J/K².mol"""
    beta = nonlinear_fit(a, b)[0]*1e-3  # conversion en J/K⁴.mol
    pi4 = np.pi**4
    theta_D = (N*k*pi4*12)/(5*beta)  # en K³
    return np.cbrt(theta_D)


def main():
    plot_fit(0, 400, cnt.squared_temperature, cnt.C_div_T, cnt.err_C_divT)


if __name__ == "__main__":
    main()
