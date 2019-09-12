
import numpy as np


def least_squares_BCES(Y1, Y2, V11, V22, V12=0, origin=False):
    """
    Make a least-squares fit for non-NaN values taking into account the errors in both rho and J variables. This implementation is based on Akritas1996 article. It is a generalization of the least-squares method. The variance of the slope is also calculated. The intersect is checked to be 0, otherwise a warning is issued.

    The fit is performed for the model
    X2i = alpha + beta * X1i + ei
    Yki = Xki + eki
    alpha = 0
    so the slope is for X2(X1) function and not the inverse.

    If origin == True, no intersect assumed. This doesn't change the lest-squares slope, but changes it's error estimate.

    Input:
        vectors of data points and errors corresponding to different embryos and ncs.

    Output:
        (beta, beta_V, alpha, alpha_V)

    """
    # Find and drop nans
    inds_not_nan = list(set(np.flatnonzero(~np.isnan(Y1))) & set(
        np.flatnonzero(~np.isnan(Y2))))
    Y1, Y2, V11, V22 = [v[inds_not_nan] for v in (Y1, Y2, V11, V22)]
    Y1m = Y1.mean()
    Y2m = Y2.mean()
    n = len(Y1)

    # Estimates for slope (beta) and intersect (alpha)
    beta = (
        np.sum((Y1 - Y1m) * (Y2 - Y2m) - V12) /
        np.sum((Y1 - Y1m)**2 - V11)
    )
    if not origin:
        alpha = (Y2m - beta * Y1m)
    else:
        alpha = 0

    # Error on the estimates
    ksi = ((Y1 - Y1m) * (Y2 - beta * Y1 - alpha) + beta * V11 - V12) / (Y1.var() - V11.mean())
    zeta = Y2 - beta * Y1 - Y1m * ksi

    beta_V = ksi.var() / n
    alpha_V = zeta.var() / n

    # T, _, _, _ = np.linalg.lstsq(slopes[:, np.newaxis], Ns, rcond=None)
    # print(beta, np.sqrt(beta_V), alpha, np.sqrt(alpha_V))
    # print('Finished!')
    return (beta, beta_V, alpha, alpha_V)
