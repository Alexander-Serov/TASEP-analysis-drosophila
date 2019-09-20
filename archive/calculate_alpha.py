

import numpy as np
from scipy.stats import shapiro

from constants import k, l, lV
from mixed_estimators import mixed_estimator_2


def alpha_over_k_J(JoK, JoKV):
    """
    Calculate mean and variance of alpha/k from J/k
    """
    alpha_over_k = (1
                    - JoK * (l - 1)
                    - np.sqrt(1 + JoK**2 * (l - 1) **
                              2 - 2 * JoK * (l + 1))
                    ) / 2
    partial_JoK = (1 / 2) * (1 - l - (2 * JoK * (-1 + l)**2 - 2 * (1 + l)) /
                             (2 * np.sqrt(1 + JoK ** 2 * (-1 + l)**2 - 2 * JoK * (1 + l))))
    partial_l = (1 / 2) * (-JoK - (-2 * JoK + 2 * JoK**2 * (-1 + l)) /
                           (2 * np.sqrt(1 + JoK**2 * (-1 + l)**2 - 2 * JoK * (1 + l))))
    alpha_over_kV = partial_l**2 * lV + partial_JoK**2 * JoKV
    return alpha_over_k, alpha_over_kV


def alpha_over_k_rho(rho, rhoV):
    """
    Calculate mean and variance of alpha/k from rho
    """
    alpha_over_k = rho / (l - rho * (l - 1))
    partial_l = -(((1 - rho) * rho) / (l - (-1 + l) * rho)**2)
    partial_rho = l / (l + rho - l * rho)**2

    alpha_over_kV = partial_l**2 * lV + partial_rho**2 * rhoV

    return alpha_over_k, alpha_over_kV


def calculate_alpha(analyses_in):
    """
    Calculates alpha/k based on experimental slopes and steady state N.
    For alpha estimates, we use rho and JoK estimates that are already averaged over differnt embryos of the same group.
    Each group corresponds to a different gene-construct-nc combination.
    After that we construct a mixed estimator based on the two following the method described in Lavancier and Rocher (2016).

    The variance of individual alpha estimators is calculated as a sum of squares of variances with partial derivatives:
    V = \sum_i (d alpha / d \theta_i)**2 * var(theta_i)
    """

    analyses = analyses_in.copy()
    significance_level = 0.05

    genes = set(analyses.gene)
    constructs = set(analyses.construct)
    ncs = set(analyses.index.get_level_values(1))

    # Load alphas obtained from J/k
    JoK = analyses.JoK
    JoKV = analyses.JoKV
    analyses['alpha_over_k_J'], analyses['alpha_over_k_JV'] = alpha_over_k_J(JoK, JoKV)
    analyses['alpha_J'] = analyses['alpha_over_k_J'] * k

    # Load alphas obtained from rho
    rho = analyses.rho
    rhoV = analyses.rhoV
    analyses['alpha_over_k_rho'], analyses['alpha_over_k_rhoV'] = alpha_over_k_rho(rho, rhoV)
    analyses['alpha_rho'] = analyses['alpha_over_k_rho'] * k

    # I use aoK to denote alpha/k
    def tau(aoK):
        # Invert alpha to get tau and recalculate into seconds
        tau = 1 / aoK / k * 60
        return tau

    refuted_shapiro_alpha = []
    refuted_shapiro_tau = []
    for gene in genes:
        for construct in constructs:
            indices = (analyses.gene == gene) & (
                analyses.construct == construct)
            for nc in ncs:
                # %% Mixed estimator for alpha/k
                # Drop nans
                T1 = drop_nan(analyses.loc[(indices, nc), 'alpha_over_k_rho'].values)
                T2 = drop_nan(analyses.loc[(indices, nc), 'alpha_over_k_J'].values)

                # Calculate the mixed estimator
                aoK_mixed, aoK_V, weights = mixed_estimator_2(T1=T1, T2=T2)

                # Save to analyses
                analyses.loc[(indices, nc), 'alpha_over_k_comb'] = aoK_mixed
                analyses.loc[(indices, nc), 'alpha_over_k_combV'] = aoK_V
                analyses.loc[(indices, nc),
                             'kappa'] = weights[0]
                analyses.loc[(indices, nc), 'alpha_comb'] = aoK_mixed * k
                analyses.loc[(indices, nc), 'alpha_combV'] = aoK_V * k**2

                # Perform Shapiro-Wilk's normality test to see whether alpha or tau are more appropriate for comparison with other conditions. 1st part
                # Save to analyses for later re-use
                data_mixed = weights[0] * T1 + weights[1] * T2
                if len(data_mixed) >= 3:
                    shap_W, shap_p = shapiro(data_mixed)
                    analyses.loc[(indices, nc), 'shap_alpha_W'] = shap_W
                    analyses.loc[(indices, nc), 'shap_alpha_p'] = shap_p
                    refuted_shapiro_alpha.append(shap_p < significance_level)

                # %% Mixed estimator for tau
                # Calculate, then save to analyses
                T1 = drop_nan(tau(analyses.loc[(indices, nc), 'alpha_over_k_rho'].values))
                T2 = drop_nan(tau(analyses.loc[(indices, nc), 'alpha_over_k_J'].values))
                tau_mixed, tau_V, weights = mixed_estimator_2(T1=T1, T2=T2)
                analyses.loc[(indices, nc), 'tau'] = tau_mixed
                analyses.loc[(indices, nc), 'tauV'] = tau_V
                analyses.loc[(indices, nc),
                             'kappa_tau'] = weights[0]

                # Perform Shapiro-Wilk's normality test for the mixed estimator
                data_mixed = weights[0] * T1 + weights[1] * T2
                if len(data_mixed) >= 3:
                    shap_W, shap_p = shapiro(data_mixed)
                    analyses.loc[(indices, nc), 'shap_tau_W'] = shap_W
                    analyses.loc[(indices, nc), 'shap_tau_p'] = shap_p
                    refuted_shapiro_tau.append(shap_p < significance_level)

                # Print the mixed estimator for the current gene-construct-nc combination
                print(
                    f'{gene}, {construct}, nc{nc}: tau = {tau_mixed:.1f} +- {np.sqrt(tau_V):.1f} s, alpha = {aoK_mixed*k:.1f} +- {np.sqrt(aoK_V*k**2):.1f} min^{-1}')

    print('\nShapiro-Wilk normality test results across all genes, constructs, ncs:')
    print(
        f'For tau, normality refuted in {np.sum(refuted_shapiro_tau)} cases out of {len(refuted_shapiro_tau)} at p = {significance_level}')
    print(
        f'For alpha, normality refuted in {np.sum(refuted_shapiro_alpha)} cases out of {len(refuted_shapiro_alpha)} at p = {significance_level}')

    return analyses


def drop_nan(T):
    return T[~np.isnan(T)]


if __name__ == '__main__':
    """Some testing code"""
    alpha_over_k_J(0.00949, 3.11e-8)
    alpha_over_k_rho(0.1, 0.0049)
