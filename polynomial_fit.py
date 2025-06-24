import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import constants as cnt
import schottky_analysis as sch
import tools


def fit_with_uncertainties(model, x_data, y_data, y_err, initial_params, n_simulations=500, bounds=(-np.inf, np.inf)):
    """
    Ajuste les paramètres d'un modèle aux données expérimentales et calcule leurs incertitudes par méthode Monte Carlo.

    Args:
        model (callable): Fonction théorique f(x, *params)
        x_data (array): Données indépendantes
        y_data (array): Données dépendantes mesurées
        y_err (array): Incertitudes sur les données y
        initial_params (list): Paramètres initiaux pour l'ajustement
        n_simulations (int): Nombre de simulations Monte Carlo
        bounds (tuple): Bornes pour les paramètres ([min], [max])

    Returns:
        dict: Dictionnaire contenant les résultats
    """
    # Fonction résidus pour les moindres carrés
    def residuals(params):
        return (y_data - model(x_data, *params)) / y_err

    # Ajustement sur les données originales
    fit_result = least_squares(
        residuals,
        initial_params,
        bounds=bounds,
        method='trf'  # Trust Region Reflective
    )

    if not fit_result.success:
        raise RuntimeError(
            "L'ajustement initial a échoué : " + fit_result.message)

    best_params = fit_result.x
    params_simulations = []

    # Simulations Monte Carlo
    for _ in range(n_simulations):
        # Génération de données bruitées
        y_sim = y_data + np.random.normal(0, y_err)

        # Fonction résidus pour la simulation
        def residuals_sim(params):
            return (y_sim - model(x_data, *params)) / y_err

        # Ajustement sur données simulées
        sim_result = least_squares(
            residuals_sim,
            best_params,  # On part des paramètres optimaux
            bounds=bounds,
            method='trf'
        )

        if sim_result.success:
            params_simulations.append(sim_result.x)

    # Calcul des incertitudes
    params_simulations = np.array(params_simulations)
    params_uncertainty = np.std(params_simulations, axis=0)

    # Calcul du chi² réduit
    chi2 = np.sum(residuals(best_params)**2)
    chi2_reduced = chi2 / (len(x_data) - len(initial_params))

    return {
        'best_params': best_params,
        'params_uncertainty': params_uncertainty,
        'chi2': chi2,
        'chi2_reduced': chi2_reduced,
        'params_simulations': params_simulations,
        'success_rate': len(params_simulations) / n_simulations
    }


def polynomial_model(x_carre, alpha, beta, gamma):
    a = alpha*(x_carre**2)
    b = beta*x_carre
    c = gamma
    return a + b + c


# results = fit_with_uncertainties(polynomial_model, cnt.squared_temperature,
    # cnt.C_div_T -
    # sch.schottky(cnt.squared_temperature,
    #            cnt.E_curve_fit, cnt.n_curve_fit)/cnt.temper,
    # cnt.err_C_divT, [0, 0.45, 1e-4])


def poly_fit(a, b, x, x_carre, y, err_y, E, n, deg=2):
    x_carre_bounded, y_bounded = tools.tab_interval(x_carre, y, a**2, b**2)
    x_bounded, erry_bounded = tools.tab_interval(x, err_y, a, b)
    schottky = sch.schottky(x_bounded, E, n)/x_bounded
    y_substracted = y_bounded - schottky
    fit = np.polyfit(x_carre_bounded, y_substracted, deg=deg)
    plt.figure()
    plt.plot(x_carre_bounded, y_substracted, "g.", label="Exp")
    plt.plot(x_carre_bounded, polynomial_model(
        x_carre_bounded, fit[0], fit[1], fit[2]), label="Fit")
    plt.errorbar(x_carre_bounded, y_substracted, erry_bounded, 0, "+g")
    plt.legend()
    plt.show()
    return fit


def main():
    pass


if __name__ == "__main__":
    main()
