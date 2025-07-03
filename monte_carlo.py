import numpy as np
import numpy.random as rd
import constants as cnt
import tools


def model_linear(x, beta, gamma):
    return beta*x + gamma


def chi2(x_data, y_data, f, beta, gamma):
    """Calculate the chi2
    x_data, y_data : array-like
    f : model 
    params = array-like of the parameters"""
    r = y_data - f(x_data, beta, gamma)
    return np.sum(r**2)


def chi2_one_iteration(x_data, y_data, f, bounds_beta, bounds_gamma):
    beta = rd.uniform(bounds_beta[0], bounds_beta[1])
    gamma = rd.uniform(bounds_gamma[0], bounds_gamma[1])
    return chi2(x_data, y_data, f, beta, gamma), beta, gamma


def minimize_chi2(x_data, y_data, f, N, bounds_beta, bounds_gamma):
    i = 0
    chi2_ini = 1e6
    beta_opt = 0
    gamma_opt = 0
    N = int(N)
    for i in range(0, N):
        chi2, beta, gamma = chi2_one_iteration(
            x_data, y_data, f, bounds_beta, bounds_gamma)
        if chi2 < chi2_ini:
            chi2_ini = chi2
            beta_opt = beta
            gamma_opt = gamma
        i += 1
    return chi2_ini, beta_opt, gamma_opt


def main():
    temperature_squared_bounded, hc_div_temp_bounded = tools.tab_interval(
        cnt.squared_temperature_HPHT, cnt.hc_div_temp_HPHT, 25, 150)
    print(minimize_chi2(temperature_squared_bounded,
          hc_div_temp_bounded, model_linear, 1e6, [0.4, 0.5], [2, 10]))


if __name__ == "__main__":
    main()
