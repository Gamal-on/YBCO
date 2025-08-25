import numpy as np
import numpy.random as rd
import algo.tools
import matplotlib.pyplot as plt


def chi2(x_data, y_data, f, params):
    """Calculate the chi2
    x_data, y_data : array-like
    f : model 
    params = array-like of the parameters"""
    r = y_data - f(x_data, *params)
    return np.sum(r**2)


def chi2_one_iteration(x_data, y_data, f, bounds, size_params):
    params = np.ones(size_params)
    for i in range(0, size_params):
        params[i] = rd.uniform(*bounds[i])
    return chi2(x_data, y_data, f, params), params


def minimize_chi2(x_data, y_data, f, N, bounds, size_params):
    """Optimisation using the minization of chi2
    x_data, y_data : arrays-like
    f : model
    N : number of iterations, integer
    bounds : tuple of array-like
    size_params : integers, number of parameters"""
    i = 0
    chi2_ini = 1e6
    params_opt = np.ones(size_params)
    N = int(N)
    for i in range(0, N):
        chi2, params = chi2_one_iteration(
            x_data, y_data, f, bounds, size_params)
        if chi2 < chi2_ini:
            chi2_ini = chi2
            params_opt = params
        i += 1
    return chi2_ini, params_opt


def monte_carlo_fitting(x_data, y_data, f, a, b, N, bounds, size_params):
    x, y = tools.tab_interval(x_data, y_data, a, b)
    chi2, params = minimize_chi2(x, y, f, N, bounds, size_params)

    # Plotting the results
    plt.figure()
    plt.plot(x, y, ".g", label="exp")
    plt.plot(x, f(x, *params), "-c", label="fit")
    plt.grid(True)
    plt.legend()
    plt.xlabel("T² (K²)")
    plt.ylabel("Hc/T (K)")
    plt.show()

    return params, chi2


def main():
    pass


if __name__ == "__main__":
    main()
