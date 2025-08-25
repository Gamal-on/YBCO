import numpy as np
import data.constants as cnt
import scipy.integrate as sci
import algo.schottky_analysis as sch

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


def model_integral_chat(x, theta, gamma):
    """Debye model with scipy.integrate.quad for numerical integration
    x : array-like, squared temperature (T^2)
    theta : float, Debye temperature
    gamma : float, constant term

    Returns the heat capacity per mole (J/mol/K) as a function of T^2.
    """
    # y = theta / sqrt(x) = theta / T
    y = theta / np.sqrt(x)
    const = 9 * cnt.N * cnt.k * 1e3  # 9R, with R = N_A * k_B, factor 1e3 for mJ

    def integrand(t):
        # More stable than (np.exp(t)-1)
        return np.exp(t) * (t**4) / (np.expm1(t)**2)

    # Integrate from 0 to y for each y
    I = np.array([sci.quad(integrand, 0, yi, epsabs=1e-8)[0]
                  for yi in np.atleast_1d(y)])
    return gamma + const * (x/(theta)**3) * I


def model_integral_schottky(x, theta, gamma, E, n):
    return model_integral_chat(x, theta, gamma) + sch.schottky(np.sqrt(x), E, n) / np.sqrt(x)


def model_einstein_schottky(x, w1, w2, F, gamma, E, n):
    """Einstein model with Schottky contribution
    x : array-like, squared temperature (T^2)
    w1 : float, Einstein frequency 1 (in K)
    w2 : float, Einstein frequency 2 (in K)
    F : float, proportion of the first Einstein mode    """
    y = np.sqrt(x)
    schottky = sch.schottky(y, E, n)/y
    mode_phonon1 = F*(w1/y)**2 * np.exp(w1/y) / (np.expm1(w1/y)**2)
    mode_phonon2 = (1-F)*(w2/y)**2 * np.exp(w2/y) / (np.expm1(w2/y)**2)
    return gamma + schottky + (mode_phonon1 + mode_phonon2) * 13 * 3 * (cnt.r/y) * 1e3


def model_einstein(x, w1, w2, F, gamma):
    """Einstein model without Schottky contribution
    x : array-like, squared temperature (T^2)
    w1 : float, Einstein frequency 1 (in K)
    w2 : float, Einstein frequency 2 (in K)
    F : float, proportion of the first Einstein mode
    gamma : float, constant term"""
    y = np.sqrt(x)
    mode_phonon1 = F*(w1/y)**2 * np.exp(w1/y) / (np.expm1(w1/y)**2)
    mode_phonon2 = (1-F)*(w2/y)**2 * np.exp(w2/y) / (np.expm1(w2/y)**2)
    return gamma + (mode_phonon1 + mode_phonon2) * 13 * 3 * (cnt.r/y) * 1e3


def main():
    pass


if __name__ == "__main__":
    main()
