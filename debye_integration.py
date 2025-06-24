import numpy as np
import matplotlib.pyplot as plt
import tools
import constants as cnt
from scipy.integrate import quad


# Step 1 : defining the function to integrate

def fonction(x):
    num = np.exp(x)*(x**4)
    denom = (np.exp(x) - 1)**2
    return num/denom


def fonction_debye(temperature, temp_debye):
    """Return the Debye function to be integrated"""
    alpha = 1/temperature
    return alpha*fonction(alpha*temp_debye)


# Step 2 : integrate the function between bounds, regarding the temperature


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

# Step 3 : Calculate the theorical beta and test it


def tab_beta(temperature):
    """Return the beta value, function of the temperature, for several possible Debye temperatures"""
    debye_temperature, integral = debye_integral(temperature)
    tab_beta = integral/(debye_temperature**3)
    return debye_temperature, 9*cnt.N*cnt.k*tab_beta


def main():
    print(tab_beta(10)[1])


if __name__ == "__main__":
    main()
