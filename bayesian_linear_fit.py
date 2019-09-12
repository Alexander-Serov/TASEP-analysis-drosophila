"""
Perform a Bayesian linear fit.
If c = False, the line must pass through the origin
See D'Agostini2005 for further details of the method.

If the fit converges to sigmaV = 0, that means that the supplied variance of the data can already explain all the variance observed. No need to attach variance to the equaiton.

Input:
prior(m, sigmaV) or prior(m, c, sigmaV) - by default is a uniform prior on all vars. Does not need to be proper

Return:
m, mV - line slope and its variance
c, cV - intersect and its variance
sigmaV, sigmaVV - equation scatter and the variance of its estimator
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize as opt

atol = 1e-16


def bayesian_linear_fit(x, y, Vx, Vy, c=True, prior=None):
    plot = 0
    sigmaV_guess = 0
    m_guess = 1
    # print(x, y)

    if c:
        phi = posterior(x, y, Vx, Vy, c=True, prior=prior)
        guess = (m_guess, 0, sigmaV_guess)  # (m, c, sigmaV)
        min = opt.minimize(phi, guess)  # , bounds=((None, None), (None, None), (0, None)))
        m_est, c_est, sigmaV_est = min.x
        mV_est, cV_est, sigmaVV_est = min.hess_inv.diagonal()

    else:
        phi = posterior(x, y, Vx, Vy, c=False, prior=prior)
        guess = (m_guess, sigmaV_guess)  # (m, sigmaV)
        min = opt.minimize(phi, guess)  # , bounds=((None, None), (0, None)))
        m_est, sigmaV_est = min.x
        mV_est, sigmaVV_est = min.hess_inv.diagonal()

        c_est = 0
        cV_est = 0
    # Find the minimum of phi (max. likelihood)

    # print(min)
    # print(np.sqrt(min.hess_inv.todense()))

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
    # print(estimates)
    # if c:
    #     res['c', 'cV'] = [c_est, cV_est]

    return estimates


def posterior(x, y, Vx, Vy, c=True, prior=None):
    """
    Returns phi = -ln(posterior), non-normalized as a function of the control parameters (m, c, Vv)
    Using a flat prior for sigmaV and intersect, and a uniform angle prior for the slope.

    posterior = p(m, c, sigmaV | x, y, Vx, Vy)
    """

    # Identify and drop nans
    inds_not_nan = list(set(np.flatnonzero(~np.isnan(x * y))))
    x, y, Vx, Vy = [ar[inds_not_nan] for ar in [x, y, Vx, Vy]]
    # # Concatenate the prior with the data set
    # x, y, Vx, Vy = [np.append(ar[inds_not_nan], el)
    #                 for ar, el in zip([x, y, Vx, Vy], [xp, yp, Vpx, Vpy])]

    # Default uniform prior
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
