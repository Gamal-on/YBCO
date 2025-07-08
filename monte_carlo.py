import numpy as np
import numpy.random as rd
import constants as cnt
import tools
import scipy.integrate as sci


def model_linear(x, beta, gamma):
    return beta*x + gamma


def model_integral(x, theta, gamma):
    """Debye model without approximating the integral
    x : array-like, squared temperature
    theta : float, Debye temperature
    gamma : float, constant term"""

    y = theta/np.sqrt(x)

    # Define the integrand function

    def integrand(t):
        """Integrand function for the integral"""
        num = np.exp(t) * (t**4)
        denom = (np.exp(t) - 1)**2
        return num / denom

    # Term integral
    integral = (1/(y**3))*sci.quad(integrand, 0, y)[0]

    return gamma + 9*cnt.N * cnt.k*integral

# Modèle avec quad (chat)

def model_integral_chat(x, theta, gamma, N):
    y = theta/np.sqrt(x)
    const = 9e3*N*cnt.k # problème unité

    def integrand(t):
        return np.exp(t)*(t**4)/(np.expm1(t)**2)  # expm1(t)=exp(t)-1, plus stable

    # on intègre de 0→y pour chaque valeur de y
    I = np.array([sci.quad(integrand, 0, yi, epsabs=1e-8)[0] for yi in np.atleast_1d(y)])
    return gamma + const * (1/(y**3))*I



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
    x_interval, y_interval = tools.tab_interval(x_data, y_data, a, b)
    return minimize_chi2(x_interval, y_interval, f, N, bounds, size_params)


def main():
    print(monte_carlo_fitting(cnt.squared_temperature_HPHT, cnt.hc_div_temp_HPHT,
          model_integral_chat, 36, 400, 5e3, ([350, 550], [0, 10], [8e23, 1e26]), 3))


if __name__ == "__main__":
    main()
