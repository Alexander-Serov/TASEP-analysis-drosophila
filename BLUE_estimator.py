import numpy as np
from numpy.linalg import inv
from tqdm import trange

from constants import l


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
    print(S, Sinv, e, T)

    # Calculate the composite estimator
    weights = e.T @ Sinv / (e.T @ Sinv @ e)
    theta_wave = weights @ T

    # kappa_opt = (x2V - cov) / (x1V + x2V - 2 * cov)
    # x = x1 * kappa_opt + x2 * (1 - kappa_opt)
    # xV = kappa_opt**2 * x1V + (1 - kappa_opt)**2 * x2V + 2 * kappa_opt * (1 - kappa_opt) * cov
    # print('BLUE', x1, x1V, x2, x2V, cov, kappa_opt, x, xV)

    # return [x, xV, kappa_opt]
    return theta_wave, weights


def mixed_estimator(T1, T2, verbose=False):
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
