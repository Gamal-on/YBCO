import numpy as np
import numpy.random as rd
import constants as cnt
import tools


def model_linear(x, beta, gamma):
    return beta*x + gamma


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


def ajustement_mc(x_data, y_data, f, a, b, N, bounds, size_params):
    x_interval, y_interval = tools.tab_interval(x_data, y_data, a, b)
    return minimize_chi2(x_interval, y_interval, f, N, bounds, size_params)


def main():
    print(ajustement_mc(cnt.squared_temperature_HPHT, cnt.hc_div_temp_HPHT,
          model_linear, 36, 144, 1e5, ([0.4, 0.5], [0, 10]), 2))


if __name__ == "__main__":
    main()
