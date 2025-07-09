import numpy as np
import constants as cnt
import matplotlib.pyplot as plt
import tools
import scipy.optimize as sc
import monte_carlo as mc
import model

# Fitting the experimental data


def curve_fit(x_data, y_data, a, b, model, bounds):
    x, y, = tools.tab_interval(x_data, y_data, a, b)
    res = sc.curve_fit(model, x, y, bounds=bounds)
    return res[0]

# Fit a curve using the BFGS method


def fit_curve_lbfgsb(x_data, y_data, a, b, model, bounds):
    """
    Ajuste une courbe expérimentale à un modèle numérique à plusieurs paramètres
    en utilisant la méthode L-BFGS-B (avec bornes).
    """
    x, y = tools.tab_interval(x_data, y_data, a, b)

    # Générer une estimation initiale au centre des bornes
    initial_guess = [(low + high) / 2 for (low, high) in bounds]

    def objective(params):
        return np.sum((y - model(x, *params))**2)

    result = sc.minimize(objective, initial_guess,
                         method='L-BFGS-B', bounds=bounds)
    return result.x, result.fun

# With ADAM


def fit_curve_adam(x_data, y_data, a, b, model, bounds, lr=0.01, epochs=1000, beta1=0.9, beta2=0.999, eps=1e-8):
    """
    Optimise une courbe expérimentale à un modèle numérique à plusieurs paramètres
    en utilisant l'algorithme Adam.
    """
    x, y = tools.tab_interval(x_data, y_data, a, b)
    params = np.array(
        [(low + high) / 2 for (low, high) in bounds], dtype=float)
    m = np.zeros_like(params)
    v = np.zeros_like(params)

    def objective(params):
        return np.sum((y - model(x, *params))**2)

    def grad(params):
        g = np.zeros_like(params)
        h = 1e-5
        for i in range(len(params)):
            params1 = params.copy()
            params2 = params.copy()
            params1[i] += h
            params2[i] -= h
            g[i] = (objective(params1) - objective(params2)) / (2 * h)
        return g

    for t in range(1, epochs + 1):
        g = grad(params)
        m = beta1 * m + (1 - beta1) * g
        v = beta2 * v + (1 - beta2) * (g ** 2)
        m_hat = m / (1 - beta1 ** t)
        v_hat = v / (1 - beta2 ** t)
        params -= lr * m_hat / (np.sqrt(v_hat) + eps)
        # Respect bounds
        for i, (low, high) in enumerate(bounds):
            params[i] = np.clip(params[i], low, high)
    return params, objective(params)


# Plotting the results


def plot_fit_gradient(x_data, y_data, a, b, model, method, bounds):
    x, y, = tools.tab_interval(x_data, y_data, a, b)
    params, chi2 = method(x_data, y_data, a, b, model, bounds)

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
