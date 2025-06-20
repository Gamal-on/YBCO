import numpy as np
import matplotlib.pyplot as plt
import tools
import fitutils as ft
import schottky_analysis
import constants

# Substraction of Schottky contribution


def schottky_substraction(E, n, x, y):
    """Return the data C/T with Schottky contribution Cs(T) substracted :
    E : scalar, n : scalar, x = temperature, y = C/T"""
    substracted = y - schottky_analysis.schottky(x, E, n)/x
    return substracted

# Linear fit with Monte Carlo


def linear_fit(a, b, x_carre, y, err_x=constants.err_squared_temperature, err_y=constants.err_sample_HC):
    """Perform a linear fit using Monte Carlo method, between a and b (bounds)
    Return the fitted values"""
    x_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    fit = ft.linfitxy(x_interval, y_interval, err_x, err_y,
                      plot=True, markercolor="g", linecolor="c")
    plt.show()
    return fit


# Debye temperature and gamma

def debye_temperature(a, b, E, n, N=78e23):
    """
    Calculate the Debye temperature and gamma from the linear fit parameters
    Returns the Debye temperature in K, gamma in J/K².mol and their respectiv errors"""
    beta, gamma, u_beta, u_gamma = linear_fit(
        a, b, E, n)*1e-3  # conversion en J
    pi4 = np.pi**4
    theta_D = (N*constants.k*pi4*12)/(5*beta)  # en K³
    u_theta_D = np.cbrt(theta_D) * u_beta/(3*beta)
    temp_debye = np.cbrt(theta_D)
    return temp_debye, gamma, u_theta_D, u_gamma


# Main function

def final(a, b, E, n, x, x_carre, y):
    """Warning : bounds are squared bounds"""
    y_substracted = schottky_substraction(E, n, x, y)
    fit = linear_fit(a, b, x_carre, y_substracted)
    return fit


def main():
    final(0, 100, constants.E_curve_fit, constants.n_curve_fit,
          constants.temperature, constants.squared_temperature, constants.C_div_T)


if __name__ == "__main__":
    main()
