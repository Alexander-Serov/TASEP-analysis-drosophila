"""
This file contains general statistical or non-statistical procedures not depending on the structure of the data
"""

import os
import shutil

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from scipy import optimize as opt
from tqdm import trange

from constants import k, kV, l, lV, tau_rev, tau_revV


def reinit_folder(folders):
    """
    Clear folder contents or create the folder if necessary.
    """
    # If string, convert to list
    if isinstance(folders, str):
        folders = [folders]

    for folder in folders:
        if not os.path.isdir(folder):
            os.makedirs(folder)
        else:
            for file in os.listdir(folder):
                file_path = os.path.join(folder, file)
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)

        print(f"Folder '{folder}' cleaned/created successfully!")


"""
This file contains data estimation and statistical calculation functions specifically written for the transcription data set.
"""


def bayesian_linear_fit(x, y, Vx, Vy, c=True, prior=None):
    """
    Perform a Bayesian linear fit for a heteroscedastic set of points with uncertainties along both axes.
    See the accompanying article and D'Agostini2005 for further details of the method.

    The method allows not only to attach uncertainties to the data points, but also to the functional relation itself.
    Basically, this means that the relation is a (normal) distribution around its mean rather than a deterministic function relating x and y.
    If the fit converges to sigmaV = 0, this means that the supplied variance of the data can already explain all the variance observed. No need to attach variance to the equaiton.

    By default, we use a flat prior for sigmaV and intersect, and a uniform angle prior for the slope.
    A custom prior function can be supplied as a parameter

    Parameters:
    x, y    -   float, array-like, lists of coordinates of 2D data points
    Vx, Vy  -   float, array-like, the corresponding of each data point
    c       -   boolean, whether the algorithm should look for a non-zero intersect. Set c = False to make the fit pass through the origin
    prior(m, sigmaV) or prior(m, c, sigmaV) - a prior function. Does not need to be proper

    Return:
    m, mV - line slope and its variance
    c, cV - intersect and its variance
    sigmaV, sigmaVV - equation scatter and the variance of its estimator
    """

    sigmaV_guess = 0
    m_guess = 1
    plot = False

    if c:
        # If not passing through the origin
        phi = posterior(x, y, Vx, Vy, c=True, prior=prior)
        guess = (m_guess, 0, sigmaV_guess)
        min = opt.minimize(phi, guess)
        m_est, c_est, sigmaV_est = min.x

        # Calculate the uncertainty on the estimates by calculating the Hessian inverse
        mV_est, cV_est, sigmaVV_est = min.hess_inv.diagonal()

    else:
        # If passing through the origin
        phi = posterior(x, y, Vx, Vy, c=False, prior=prior)
        guess = (m_guess, sigmaV_guess)
        min = opt.minimize(phi, guess)
        m_est, sigmaV_est = min.x
        mV_est, sigmaVV_est = min.hess_inv.diagonal()

        c_est = 0
        cV_est = 0

    if plot:
        plt.figure(clear=True, num=1)
        plt.errorbar(x, y, xerr=np.sqrt(Vx), yerr=np.sqrt(Vy), fmt='.', elinewidth=0.5)
        plt.ylim(ymin=0)
        plt.xlim(xmin=0)

        xfit = [np.nanmin(x), np.nanmax(x)]
        yfit = np.array(xfit) * m_est + c_est
        plt.plot(xfit, yfit, 'r')

    # Construct the result dictionary
    estimates = {'m': m_est, 'mV': mV_est, 'c': c_est,
                 'cV': cV_est, 'sigmaV': sigmaV_est, 'sigmaVV': sigmaVV_est}

    return estimates


def posterior(x, y, Vx, Vy, c=True, prior=None):
    """
    Returns phi = -ln(posterior), non-normalized as a function of the control parameters (m, c, Vv)
    By default, we use a flat prior for sigmaV and intersect, and a uniform angle prior for the slope.

    posterior = p(m, c, sigmaV | x, y, Vx, Vy)
    """

    # Drop nans
    inds_not_nan = list(set(np.flatnonzero(~np.isnan(x * y))))
    x, y, Vx, Vy = [ar[inds_not_nan] for ar in [x, y, Vx, Vy]]

    # Default prior
    if not prior:
        if c:
            def prior(m, c, sigmaV):
                return 1 / np.pi / (m**2 + 1)
        else:
            def prior(m, sigmaV):
                return 1 / np.pi / (m**2 + 1)

    def ln_likelihood(m, c, sigmaV):
        return -np.sum(np.log(sigmaV**2 + Vy + m**2 * Vx) / 2
                       + (y - m * x - c)**2 / 2 / (sigmaV**2 + Vy + m**2 * Vx)
                       )

    def phi(params):
        m, c, sigmaV = params
        phi = -ln_likelihood(m, c, sigmaV) - np.log(prior(m, c, sigmaV))
        return phi

    def phi_no_c(params):
        m, sigmaV = params
        phi = -ln_likelihood(m, 0, sigmaV) - np.log(prior(m, sigmaV))
        return phi

    if c:
        return phi
    else:
        return phi_no_c


def BLUE_estimator(x1, x1V, x2, x2V, cov=0):
    """
    This function is updated to implement the correlated estimator from Keller, Olkin (2004) for the case of k=2.
    Diregard the old description below.


    # A BLUE estimator for combining estimates
    # x = k* x1 + (1-k) * x2.
    # Assumes x1 and x2 have the same means.
    """

    # The sample covariance matrix
    S = np.array([[x1V, cov],
                  [cov, x2V]])
    Sinv = np.linalg.inv(S)
    e = np.ones((2, 1))
    T = np.array([[x1, x2]]).T
    # print(S, Sinv, e, T)

    # Calculate the composite estimator
    weights = e.T @ Sinv / (e.T @ Sinv @ e)
    mixed_estimator = weights @ T

    mixedV = 1 / (e.T @ Sinv @ e)[0, 0]

    # kappa_opt = (x2V - cov) / (x1V + x2V - 2 * cov)
    # x = x1 * kappa_opt + x2 * (1 - kappa_opt)
    # xV = kappa_opt**2 * x1V + (1 - kappa_opt)**2 * x2V + 2 * kappa_opt * (1 - kappa_opt) * cov
    # print('BLUE', x1, x1V, x2, x2V, cov, kappa_opt, x, xV)

    # return [x, xV, kappa_opt]
    return mixed_estimator[0, 0], mixedV, weights[0, 0]


def mixed_estimator_2(T1, T2, verbose=False):
    """
    Based on the Lavancier and Rochet (2016) article.

    The method combines two series of estimates of the same quantity taking into account their correlations. The individual measureements are assumed independent. The current implementation works only for point estimates.

    The main result corresponds to Eq. (11) from the article.
    Its variance is the equation after Eq. (9). [equation not checked]
    """

    B = 1000  # bootstrap repetitions

    # Drop nans
    not_nans = np.logical_or(np.isnan(T1), np.isnan(T2))
    T1, T2 = T1[~not_nans], T2[~not_nans]

    n = len(T1)
    # Return nan if no samples
    # If one sample, return simple average with no variance
    if n == 0:
        return np.nan, np.nan, np.array([np.nan, np.nan])
    elif n == 1:
        # print(T1)
        return T1[0] / 2 + T2[0] / 2, np.nan, np.array([0.5, 0.5])

    # Calculate the estimators for the data set. This is the input data for the rest
    T1_data_median = np.median(T1)
    T2_data_median = np.median(T2)

    # Estimate the covariance sigma matrix with bootstrap (with replacement, as described in the article)
    sigma = np.zeros((2, 2))
    for b in range(B):
        T1_sample = np.random.choice(T1, size=n, replace=True)
        T2_sample = np.random.choice(T2, size=n, replace=True)
        # print('T1', T1_sample)
        T1_sample_median = np.median(T1_sample)
        T2_sample_median = np.median(T2_sample)
        sigma += np.array([[(T1_sample_median - T1_data_median)**2,
                            (T1_sample_median - T1_data_median) * (T2_sample_median - T2_data_median)],
                           [(T1_sample_median - T1_data_median) * (T2_sample_median - T2_data_median),
                            (T2_sample_median - T2_data_median)**2]])
    sigma /= B

    # print(n, sigma)

    # Calculate the mixed estimator
    I = np.array([[1, 1]]).T
    T = np.array([[T1_data_median, T2_data_median]]).T
    weights = inv(I.T @ inv(sigma) @ I) @ I.T @ inv(sigma)
    mixed_estimator = (weights @ T)[0, 0]
    mixedV = (inv(I.T @ inv(sigma) @ I))[0, 0]

    if verbose:
        print('weights', weights)
        print(mixed_estimator, '+-', np.sqrt(mixedV))

    return mixed_estimator, mixedV, np.squeeze(weights)


def mixed_estimator_4(T1, T2, T3=None, T4=None, verbose=False):
    """
    Based on the Lavancier and Rochet (2016) article.

    The method combines four series of estimates of the same quantity taking into account their correlations. The individual measureements are assumed independent. The current implementation works only for point estimates.

    The main result corresponds to Eq. (11) from the article.
    Its variance is the equation after Eq. (9). [equation not checked]
    """

    B = 1000  # bootstrap repetitions

    # Drop nans
    Ts = [T1, T2, T3, T4]
    not_nans = np.isnan(T1)
    for t in [T2, T3, T4]:
        not_nans = np.logical_or(not_nans, np.isnan(t))

    for t in [T1, T2, T3, T4]:
        t = t[~not_nans]    # will work because they are aliases, not copies

    n = len(T1)

    # Calculate the estimators for the data set. This is the input data for the rest
    T_data_medians = np.full((4, 1), np.nan)
    for i, t in enumerate(Ts):
        T_data_medians[i] = np.median(t)

    I = np.ones((4, 1))
    T = T_data_medians  # np.array([[T1_data_median, T2_data_median]]).T
    # not_nans = np.logical_or(np.isnan(T1), np.isnan(T2))
    # T1, T2 = T1[~not_nans], T2[~not_nans]

    # Return nan if no samples
    # If one sample, return simple average with no variance
    if n == 0:
        return np.nan, np.nan, np.full((1, 4), np.nan)
    elif n == 1:
        weights = np.full((1, 4), 1 / 4)
        return weights @ T, np.nan, weights

    # Estimate the covariance sigma matrix with bootstrap (with replacement, as described in the article)
    sigma = np.zeros((4, 4))
    T_sample_medians = np.full((4, 1), np.nan)
    for b in range(B):
        for i, t in enumerate(Ts):
            # T_sample[i] = np.random.choice(t, size=n, replace=True)
            T_sample_medians[i] = np.median(np.random.choice(t, size=n, replace=True))

        # T1_sample = np.random.choice(T1, size=n, replace=True)
        # T2_sample = np.random.choice(T2, size=n, replace=True)
        # print('T1', T1_sample)
        # T1_sample_median = np.median(T1_sample)
        # T2_sample_median = np.median(T2_sample)

        # столбец на строку
        sigma += (T_sample_medians - T_data_medians) @ (T_sample_medians - T_data_medians).T

        # sigma += np.array([[(T1_sample_median - T1_data_median)**2,
        #                     (T1_sample_median - T1_data_median) * (T2_sample_median - T2_data_median)],
        #                    [(T1_sample_median - T1_data_median) * (T2_sample_median - T2_data_median),
        #                     (T2_sample_median - T2_data_median)**2]])
    sigma /= B

    # print(n, sigma)

    # Calculate the mixed estimator

    weights = inv(I.T @ inv(sigma) @ I) @ I.T @ inv(sigma)
    mixed_estimator = (weights @ T)[0, 0]
    mixedV = (inv(I.T @ inv(sigma) @ I))[0, 0]

    if verbose:
        print('weights', weights)
        print(mixed_estimator, '+-', np.sqrt(mixedV))

    return mixed_estimator, mixedV, np.squeeze(weights)


def welchs_test(x1Mean, x1V, n1, x2Mean, x2V, n2):
    """
    Can accept arrays as input
    """
    # Test statistics
    t = (x1Mean - x2Mean) / \
        np.sqrt(x1V / n1 + x2V / n2)
    nu = ((x1V / n1 + x2V / n2)**2
          / (x1V**2 / n1 / (n1 - 1) + x2V**2 / n2 / (n2 - 1))
          )
    # p-value for additivity
    if not np.isnan(t):
        p = 2 * (student.cdf(-np.abs(t), df=nu))
    else:
        p = np.nan

    return p, t, nu


def welchs_test_arrays(x1, x2):
    """
    A spearate implementation that automatically calculates means, variances and the number of samples from the input lists of values.

    Not currently used in the code because the number of samples is calculated earlier, when the estimators were combined.
    It is kept as a separate number for the mixed estimators.
    """

    x1, x2 = x1[~np.isnan(x1)], x2[~np.isnan(x2)]
    n1, n2 = len(x1), len(x2)

    x1Mean = np.mean(x1)
    x2Mean = np.mean(x2)

    x1V = np.var(x1, ddof=1)
    x2V = np.var(x2, ddof=1)

    return welchs_test(x1Mean, x1V, n1, x2Mean, x2V, n2)


def set_figure_size(num, rows, page_width_frac, height_factor=1.0, clear=True):
    pagewidth_in = 6.85
    font_size = 8
    dpi = 100

    figsize = np.asarray([1.0,  rows *
                          height_factor]) * page_width_frac * pagewidth_in  # in inches

    # Set default font size
    matplotlib.rcParams.update({'font.size': font_size})

    # Enable LaTeX and set font to Helvetica
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}',    # helvetica font
        r'\usepackage{sansmath}',   # math-font matching  helvetica
        r'\sansmath'                # actually tell tex to use it!
        r'\usepackage{siunitx}',    # micro symbols
        r'\sisetup{detect-all}',    # force siunitx to use the fonts
    ]

    # Create and return figure handle
    fig = plt.figure(num, clear=clear)
    fig.clear()

    # Set figure size and dpi
    fig.set_dpi(dpi)
    fig.set_figwidth(figsize[0])
    fig.set_figheight(figsize[1])

    # fig.tight_layout()

    # Return figure handle
    return (fig)


# %% Theoretical TASEP estimates
alpha_over_k_MC = 1 / (1 + np.sqrt(l))
alpha_over_k_MCV = alpha_over_k_MC**2 / 2 / np.sqrt(l)
alpha_over_k_abortive = 1 / k / tau_rev
alpha_over_k_abortiveV = alpha_over_k_abortive**2 * (kV / k**2 + tau_revV / tau_rev**2)

f'alpha_MC/k = {alpha_over_k_MC:.3f} +- {np.sqrt(alpha_over_k_MCV):.3f}'
f'alpha_abrt/k = {alpha_over_k_abortive:.4f} +- {np.sqrt(alpha_over_k_abortiveV):.4f}'


def rho_LD(alphaT):
    return alphaT * l / (1 + alphaT * (l - 1))


def rho_HD(betaT):
    return (1 - betaT)


def rho_MC():
    sq = np.sqrt(l)
    rho = sq / (1 + sq)

    partial_l = 1 / 2 / sq / (1 + sq)**2
    rhoV = partial_l**2 * lV
    # print(f'rho_MC: {rho} +- {np.sqrt(rhoV)}')
    return [rho, rhoV]


def J_over_k_LD(alphaT):
    return alphaT * (1 - alphaT) / (1 + alphaT * (l - 1))


def J_over_k_HD(betaT):
    return betaT * (1 - betaT) / (1 + betaT * (l - 1))


def J_over_k_MC():
    sq = np.sqrt(l)
    JoK = (1 + sq)**(-2)
    partial_l = - (1 + sq)**(-3) / sq
    JoKV = partial_l**2 * lV
    # print(f'J/k_MC: {JoK} +- {np.sqrt(JoKV)}')
    return [JoK, JoKV]


if __name__ == "__main__":
    rho_MC()
    J_over_k_MC()

    rho_LD(alpha_over_k_MC)

    J_over_k_LD(alpha_over_k_MC)

    3.9e-3 / J_over_k_LD(alpha_over_k_MC)
