import numpy as np
import numpy.random as rd
import constants as cnt
import tools


def model_linear(x, beta, gamma):
    return beta * x + gamma


def chi2(y_data, y_model):
    """Calcule le chi² sans pondération (incertitude uniforme)"""
    return np.sum((y_data - y_model)**2)


def minimize_chi2(x_data, y_data, model, N, bounds):
    """
    Parameters:
    x_data, y_data : data
    model : function
    N : number of iterations
    bounds : paramters intervals [(min1, max1), (min2, max2), ...]
    """
    # Convertir les bornes en format numpy
    bounds = np.array(bounds)
    size_params = bounds.shape[0]

    # Génération de tous les paramètres aléatoires
    params_matrix = np.column_stack([
        rd.uniform(low=low, high=high, size=int(N))
        for (low, high) in bounds])

    # Calcul des prédictions du modèle pour tous les jeux de paramètres
    y_models = np.array([
        model(x_data, *params) for params in params_matrix])

    # Calcul du chi² pour toutes les simulations (vectorisé)
    chi2_values = np.sum((y_models - y_data)**2, axis=1)

    # Trouver le minimum
    min_idx = np.argmin(chi2_values)
    return chi2_values[min_idx], params_matrix[min_idx]


def ajustement_mc(x_data, y_data, model, a, b, N, bounds):
    """Sélection de l'intervalle et optimisation"""
    x_interval = x_data[a:b+1]
    y_interval = y_data[a:b+1]
    return minimize_chi2(x_interval, y_interval, model, N, bounds)


def main():
    print(ajustement_mc(cnt.squared_temperature_HPHT,
          cnt.hc_div_temp_HPHT, model_linear, 36, 144, 1e7, [(0, 1), (0, 10)]))


if __name__ == "__main__":
    main()
