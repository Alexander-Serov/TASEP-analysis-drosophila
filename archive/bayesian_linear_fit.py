

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize as opt


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
