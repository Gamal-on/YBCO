import numpy as np
import matplotlib.pyplot as plt
import tools
import constants as cnt
from scipy.integrate import quad


def fonction(x):
    num = np.exp(x)*(x**4)
    denom = (np.exp(x) - 1)**2
    return num/denom


def fonction_debye(temperature, temp_debye):
    """Return the Debye function to be integrated"""
    alpha = 1/temperature
    return alpha*fonction(alpha*temp_debye)


def integrate_between_bounds(x, y, a, b):
    # Select indices where x is within [a, b]
    mask = (x >= a) & (x <= b)
    x_selected = x[mask]
    y_selected = y[mask]
    # Use the trapezoidal rule for integration
    return np.trapz(y_selected, x_selected)


def progressive_integration(x, y):
    """Return a list of the integrated values y vs x"""
    integration = []
    for val in x:
        aire = integrate_between_bounds(x, y, 0, val)
        integration.append(aire)
    return integration


def debye_integral(temperature):
    """Return a tuple of array-like objet : the chosen interval with the debye temperature and the integration"""
    debye_temperature = np.arange(0.1, 700, 0.1)
    y = fonction_debye(temperature, debye_temperature)
    integral = progressive_integration(debye_temperature, y)
    return debye_temperature, integral


def tab_beta(temperature):
    debye_temprature, integral = debye_integral(temperature)
    tab = []
    for x, y in debye_temprature, integral:
        tab.append(9*cnt.N*cnt.k * y/(x**3))
    return x, tab


def integral(y):
    if y <= 0:
        return 0

    def integrand(x): return (x**4 * np.exp(x)) / (np.exp(x) - 1)**2
    result, _ = quad(integrand, 0, y)
    return result


def F(theta_D, T_i, beta, factor):
    factor = (12/5)*np.pi**4
    y = theta_D / T_i
    return beta * theta_D**3 - factor * I(y)
