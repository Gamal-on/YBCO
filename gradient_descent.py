import numpy as np
import constants as cnt
import matplotlib.pyplot as plt
import tools
import scipy.optimize as sc
import monte_carlo as mc


# Fitting the experimental data

def curve_fit(x_data, y_data, a, b, model, bounds):
    x, y, = tools.tab_interval(x_data, y_data, a, b)
    res = sc.curve_fit(model, x, y, bounds=bounds)
    return res[0]

# Plotting the results


def plot_fit_gradient(x_data, y_data, a, b, model, bounds):
    x, y, = tools.tab_interval(x_data, y_data, a, b)
    params = curve_fit(x_data, y_data, a, b, model, bounds)

    # Plot
    plt.figure()
    plt.plot(x, y, ".g", label="exp")
    plt.plot(x, model(x, *params), "-c", label="fit")
    plt.grid(True)
    plt.legend()
    plt.xlabel("T² (K²)")
    plt.ylabel("Hc/T (K)")
    plt.show()

    # Print Monte Carlo chi2
    chi2 = mc.chi2(x, y, model, params)

    return params, chi2


def main():
    pass


if __name__ == "__main__":
    main()
