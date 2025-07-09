import numpy as np
import constants as cnt
import tools
import scipy.integrate as sci
import schottky_analysis as sch

# Al the models are defined like y(x) = C/T(T²)

# Linear model


def model_linear(x, beta, gamma):
    """Linear model for fitting"""
    return beta * x + gamma

# Non linear beta model, with the Schottky contribution (acoustic)


def model_beta_schottky(x, beta, gamma, E, n):
    """Non linear model for fitting, with Schottky contribution
    x : array-like, squared temperature
    beta : float, linear coefficient
    gamma : float, constant term
    E : float, Schottky energy
    n : float, Schottky proportionality factor"""
    cs = sch.schottky(np.sqrt(x), E, n)
    return gamma + beta * x + cs / np.sqrt(x)

# Non linear model, with the Schottky contribution (optic)


def model_alpha_schottky(x, beta, gamma, E, n, alpha):
    """Non linear model for fitting, with Schottky contribution
    x : array-like, squared temperature
    beta : float, linear coefficient
    gamma : float, constant term
    E : float, Schottky energy
    n : float, Schottky proportionality factor
    alpha : float, optic coefficient
    Return y(x) = c(x) + beta x + gamma + alpha x²"""
    cs = sch.schottky(np.sqrt(x), E, n)
    return gamma + beta * x + alpha * (x**2) + (cs / np.sqrt(x))

# Polynomial model


def model_polynomial(x, beta, gamma, alpha):
    return gamma + beta * x + alpha * (x**2)

# Model with the integral with quad


def model_integral_chat(x, theta, gamma, N):
    """Debye model with scipy.integrate.quad for numerical integration
    x : array-like, squared temperature
    theta : float, Debye temperature
    gamma : float, constant term
    N : number of atoms"""
    y = theta/np.sqrt(x)
    const = 9e3*N*cnt.k  # problème unité

    def integrand(t):
        # expm1(t)=exp(t)-1, plus stable
        return np.exp(t)*(t**4)/(np.expm1(t)**2)

    # on intègre de 0→y pour chaque valeur de y
    I = np.array([sci.quad(integrand, 0, yi, epsabs=1e-8)[0]
                 for yi in np.atleast_1d(y)])
    return gamma + const * (1/(y**3))*I


def main():
    pass


if __name__ == "__main__":
    main()
