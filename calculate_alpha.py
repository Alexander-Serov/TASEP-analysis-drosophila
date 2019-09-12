import numpy as np
import pandas as pd
from scipy.stats import shapiro

from BLUE_estimator import BLUE_estimator, mixed_estimator
from constants import L, k, kV, l, lV
from theoretical_phase_parameters import alpha_over_k_MC


def alpha_over_k_J(JoK, JoKV):
    """Inverse function for alpha/k from J/k"""
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
    """Inverse function for alpha/k from J/k"""
    alpha_over_k = rho / (l - rho * (l - 1))
    partial_l = -(((1 - rho) * rho) / (l - (-1 + l) * rho)**2)
    partial_rho = l / (l + rho - l * rho)**2

    alpha_over_kV = partial_l**2 * lV + partial_rho**2 * rhoV

    return alpha_over_k, alpha_over_kV


def calculate_alpha(analyses_in, I_est):
    """Calculates alpha/k based on experimental slopes and steady state N.
    For alpha estimates, we use rho and JoK estimates that are already averaged over differnt embryos of the same group. Each group corresponds to a different gene-construct-nc combination.
    After that constructs a minimal variance estimator based on the two.
    Requires the calibration coefficient I_est in fluo/pol

    # The variance of alpha is estimated through a sum of squares of variances with partial derivatives:
    # V = \sum (d alpha / d \theta_i)**2 * var(theta_i)
    """
    analyses = analyses_in.copy()
    sign_level = 0.05

    # From J/k
    JoK = analyses.JoK
    JoKV = analyses.JoKV
    analyses['alpha_over_k_J'], analyses['alpha_over_k_JV'] = alpha_over_k_J(JoK, JoKV)
    analyses['alpha_J'] = analyses['alpha_over_k_J'] * k

    # From rho
    rho = analyses.rho
    rhoV = analyses.rhoV
    analyses['alpha_over_k_rho'], analyses['alpha_over_k_rhoV'] = alpha_over_k_rho(rho, rhoV)
    analyses['alpha_rho'] = analyses['alpha_over_k_rho'] * k
    # print(analyses)

    # def cov_unbiased(x, y):
    #     nans = np.isnan(x * y)
    #     if np.sum(~nans) > 0:
    #         r = np.cov(x[~nans], y[~nans], ddof=1)
    #         return r[0, 1]
    #     else:
    #         return np.nan
    #
    # # Estimate average alpha by group taking into account correlations between the 2 estimates
    # grouped = analyses.groupby(by=['gene', 'construct', 'nc'], as_index=True)
    # results = pd.DataFrame()
    # # print(grouped.count()[['alpha_over_k_J', 'alpha_over_k_rho']])
    #
    # # alpha / k
    # aRmean = grouped.alpha_over_k_rho.mean()
    # aJmean = grouped.alpha_over_k_J.mean()
    # # aRV = grouped['alpha_over_k_rho'].var(ddof=1)
    # # aJV = grouped['alpha_over_k_J'].var(ddof=1)
    # # cov = grouped.apply(lambda group: cov_unbiased(
    # #     group['alpha_over_k_J'], group['alpha_over_k_rho']))
    # # print(cov)
    # # alpha_over_k_comb, alpha_over_k_combV, kappa_opt = BLUE_estimator(
    # #     aJmean, aJV, aRmean, aRV, cov)
    # # improve_k = np.min([aRV, aJV], axis=0)
    # # improve_k = (alpha_over_k_combV - improve_k) / improve_k * 100

    def tau(aoK):
        tau = 1 / aoK / k * 60
        return tau

    genes = set(analyses.gene)
    constructs = set(analyses.construct)
    ncs = set(analyses.index.get_level_values(1))
    refuted_shapiro_alpha = []
    refuted_shapiro_tau = []
    for gene in genes:
        for construct in constructs:
            for nc in ncs:
                # print(nc)
                indices = (analyses.gene == gene) & (
                    analyses.construct == construct)
    # for key, data in grouped:
        # print
        # print(key, data.rho.values)

                # Mixed estimator for alpha/k
                T1 = drop_nan(analyses.loc[(indices, nc), 'alpha_over_k_rho'].values)
                T2 = drop_nan(analyses.loc[(indices, nc), 'alpha_over_k_J'].values)
                aoK_mixed, aoK_V, weights = mixed_estimator(T1=T1, T2=T2)

                analyses.loc[(indices, nc), 'alpha_over_k_comb'] = aoK_mixed
                analyses.loc[(indices, nc), 'alpha_over_k_combV'] = aoK_V
                analyses.loc[(indices, nc),
                             'kappa'] = weights[0]

                # Perform Shapiro-Wilks normality test for the mixed estimator
                data_mixed = weights[0] * T1 + weights[1] * T2
                if len(data_mixed) >= 3:
                    # print(T1)
                    shap_W, shap_p = shapiro(data_mixed)
                    analyses.loc[(indices, nc), 'shap_alpha_W'] = shap_W
                    analyses.loc[(indices, nc), 'shap_alpha_p'] = shap_p
                    refuted_shapiro_alpha.append(shap_p < sign_level)

                    # Mixed estimator for alpha
                    # a_mixed, a_V, weights = mixed_estimator(
                    #     T1=analyses.loc[(indices, nc), 'alpha_over_k_rho'].values * k, T2=analyses.loc[(indices, nc), 'alpha_over_k_J'].values * k)
                analyses.loc[(indices, nc), 'alpha_comb'] = aoK_mixed * k
                analyses.loc[(indices, nc), 'alpha_combV'] = aoK_V * k**2
                # analyses.loc[(indices, nc),
                #              'kappa_alpha'] = weights[0]

                # Mixed estimator for alpha/k
                T1 = drop_nan(tau(analyses.loc[(indices, nc), 'alpha_over_k_rho'].values))
                T2 = drop_nan(tau(analyses.loc[(indices, nc), 'alpha_over_k_J'].values))
                tau_mixed, tau_V, weights = mixed_estimator(T1=T1, T2=T2)
                analyses.loc[(indices, nc), 'tau'] = tau_mixed
                analyses.loc[(indices, nc), 'tauV'] = tau_V
                analyses.loc[(indices, nc),
                             'kappa_tau'] = weights[0]

                # Perform Shapiro-Wilks normality test for the mixed estimator
                data_mixed = weights[0] * T1 + weights[1] * T2
                if len(data_mixed) >= 3:
                    # print(T1)
                    shap_W, shap_p = shapiro(data_mixed)
                    analyses.loc[(indices, nc), 'shap_tau_W'] = shap_W
                    analyses.loc[(indices, nc), 'shap_tau_p'] = shap_p
                    refuted_shapiro_tau.append(shap_p < sign_level)

                # print(np.sqrt(tau_V))
                # print(tau)
                # print(f'{tau: .1f}')
                print(
                    f'{gene}, {construct}, nc{nc}: tau = {tau_mixed:.1f} +- {np.sqrt(tau_V):.1f} s, alpha = {aoK_mixed*k:.1f} +- {np.sqrt(aoK_V*k**2):.1f} min^{-1}')

    print('\nShapiro-Wilk normality test results across all genes, constructs, ncs:')
    print(
        f'For tau, normality refuted in {np.sum(refuted_shapiro_tau)} cases out of {len(refuted_shapiro_tau)} at p = {sign_level}')
    print(
        f'For alpha, normality refuted in {np.sum(refuted_shapiro_alpha)} cases out of {len(refuted_shapiro_alpha)} at p = {sign_level}')

    return analyses

    # # alpha, in min^-1
    # aRmean = grouped.alpha_rho.mean()
    # aJmean = grouped.alpha_J.mean()
    # aRV = grouped.alpha_rho.var(ddof=1)
    # aJV = grouped.alpha_J.var(ddof=1)
    # cov = grouped.apply(lambda group: cov_unbiased(
    #     group.alpha_J, group.alpha_rho))
    # alpha_comb, alpha_combV, kappa_opt = BLUE_estimator(
    #     aJmean, aJV, aRmean, aRV, cov)
    # improve = np.min([aRV, aJV], axis=0)
    # improve = (alpha_combV - improve) / improve * 100
    #
    # # Estimate tau
    # tau = 1 / alpha_over_k_comb / k * 60
    # tauV = tau**2 * (kV / k**2 + alpha_over_k_combV / alpha_over_k_comb**2)
    # # print(tau, tauV)
    #
    # # Saving to analyses
    # genes = set(analyses.gene)
    # constructs = set(analyses.construct)
    # ncs = set(analyses.index.get_level_values(1))
    # for gene in genes:
    #     for construct in constructs:
    #         for nc in ncs:
    #             indices = (analyses.gene == gene) & (
    #                 analyses.construct == construct)
    #             analyses.loc[(indices, nc),
    #                          'kappa'] = kappa_opt.loc[(gene, construct, nc)]
    #             # print(kappa_opt)
    #             # print(analyses.loc[:,
    #             #                    'kappa'])
    #             # analyses.loc[(indices, nc), 'alpha_origin'] = alpha_origin.loc[(
    #             #     gene, construct, nc)]
    #             analyses.loc[(indices, nc), 'alpha_over_k_comb'] = alpha_over_k_comb.loc[(
    #                 gene, construct, nc)]
    #             analyses.loc[(indices, nc), 'alpha_over_k_combV'] = alpha_over_k_combV.loc[(
    #                 gene, construct, nc)]
    #
    #             analyses.loc[(indices, nc), 'alpha_comb'] = alpha_cur = alpha_comb.loc[(
    #                 gene, construct, nc)]
    #             analyses.loc[(indices, nc), 'alpha_combV'] = alpha_curV = alpha_combV.loc[(
    #                 gene, construct, nc)]
    #
    #             analyses.loc[(indices, nc), 'tau'] = tau_cur = tau.loc[(
    #                 gene, construct, nc)]
    #             analyses.loc[(indices, nc), 'tauV'] = tau_curV = tauV.loc[(
    #                 gene, construct, nc)]
    #
    #             # print(
    #             #     f'Combining rho and J estimates for alpha/k for {gene} {construct} nc{nc} reduced variance by {improve_k.loc[(gene, construct, nc)]:.0f}%')
    #             #
    #             # print(
    #             #     f'Combining rho and J estimates for alpha for {gene} {construct} nc{nc} reduced variance by {improve.loc[(gene, construct, nc)]:.0f}%')
    #
    #             # if nc == 13:
    #             print(
    #                 f'{gene}, {construct}, nc{nc}: tau = {tau_cur:.1f} +- {np.sqrt(tau_curV):.1f} s, alpha = {alpha_cur:.1f} +- {np.sqrt(alpha_curV):.1f} min^{-1}')
    #
    # # print(analyses)
    # # results['kappa'] = kappa_opt
    # # results['alpha_origin'] = alpha_origin
    # # results['alpha_over_k_comb'] = alpha_over_k_comb
    # # results['alpha_over_k_combV'] = alpha_over_k_combV
    # # print(np.sqrt(alpha_over_k_combV))
    #
    # # print(results)
    # # analyses = pd.merge(analyses, results, left_on=[
    # #     'gene', 'construct', 'nc'], right_on=['gene', 'construct', 'nc'])

    return analyses

    # %%
    kElong = 4e3 / 60  # bp/s
    aoK = np.array([alpha_over_k_J(JoK=0.0036), alpha_over_k_rho(rho=0.18)])
    print(aoK)
    1 / aoK / kElong


def drop_nan(T):
    return T[~np.isnan(T)]


if __name__ == '__main__':
    alpha_over_k_J(0.00949, 3.11e-8)
    alpha_over_k_rho(0.1, 0.0049)
