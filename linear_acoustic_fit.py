import numpy as np
import matplotlib.pyplot as plt
import tools
import fitutils as ft
import schottky_analysis
import constants as cnt

# Substraction of Schottky contribution


def schottky_substraction(E, n, x, y):
    """Return the data C/T with Schottky contribution Cs(T) substracted :
    E : scalar, n : scalar, x = temperature, y = C/T"""
    substracted = y - schottky_analysis.schottky(x, E, n)/x
    return substracted

# Linear fit with Monte Carlo


def linear_fit(a, b, x_carre, y, err_x, err_y):
    """Perform a linear fit using Monte Carlo method,
    a, b : bounds, squared temperature
    x_carre : array like
    y : y_data
    err_x, erry : array-like,
    Return the fitted values"""
    x_interval, y_interval = tools.tab_interval(x_carre, y, a, b)
    plt.figure()
    fit = ft.linfitxy(x_interval, y_interval, err_x[0:len(x_interval)], err_y[0:len(y_interval)],
                      plot=True, markercolor="g", linecolor="c")
    plt.xlabel("Temperature (K)")
    plt.ylabel("C/T (mJ/K²/mol)")
    plt.show()
    return fit


# Main function

def final(a, b, E, n, x, x_carre, y, err_x, err_y):
    """Substract the Schottky contribution and return the fitted paramters of C/T(T²). Warning : bounds are squared bounds"""
    y_substracted = schottky_substraction(E, n, x, y)
    fit = linear_fit(a, b, x_carre, y_substracted, err_x, err_y)
    return fit


def main():
    pass


if __name__ == "__main__":
    main()
