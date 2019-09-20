import numpy as np

from bayesian_linear_fit import bayesian_linear_fit
from constants import LV, L, k, kV


def calculate_free_travel_time(analyses_in):
    """
    Estimate the free passage time by fitting the non-scaled current-density diagram across all ncs for each gene and construct individually.

    The fitting procedure takes into account uncertainties in data points along both axes.
    the fit line must pass through 0.

    A Gaussian distribution of the slope values is used a prior.
    The width of the distribution `TV_prior` is defined by the a priori error in the elongation rate k and effective gene length L.
    Those are set in `constants.py`
    """

    analyses = analyses_in.copy()
    genes = set(analyses.gene)

    # Add columns to the analyses table
    columns = ['T', 'TV', 'L', 'LV']
    for column in columns:
        analyses[column] = np.nan

    # Make a prior on the travel time
    T_prior = L / k
    TV_prior = T_prior**2 * (LV / L**2 + kV / k**2)

    def prior_func(m, sigmaV):
        """
        The prior must be a function of 2 parameters, so sigmaV in the definition is necessary.
        But this particular prior does not depend on sigmaV
        This is a conjugated prior in the form of just the Gaussian distribution for the slope.
        """
        return np.exp(-(m - T_prior)**2 / 2 / TV_prior) * (2 * np.pi * TV_prior)**(-1 / 2)

    # Fit the diagram individually for each gene, but all constructs and ncs mixed
    # This makes sense because all constructs and nc supposedly have the same length and polymerase elongation rate
    for gene in genes:
        analyses_filtered = analyses[(analyses.gene == gene)]
        indices = analyses_filtered.index

        Ns = analyses_filtered['max'].values
        NsV = analyses_filtered['maxV'].values
        slopes = analyses_filtered['slope'].values
        slopesV = analyses_filtered['slopeV'].values

        # Fit N/s passing through zero. Give result in minutes
        fit = bayesian_linear_fit(x=slopes, Vx=slopesV, y=Ns, Vy=NsV, c=False, prior=prior_func)
        T, TV = [fit[key] for key in ['m', 'mV']]

        # Estimate L from T and k
        L_est = k * T
        LVest = T**2 * kV + k**2 * TV

        # Round up and print
        L_rnd, Lstd_rnd = np.round([L_est, np.sqrt(LVest)], decimals=-1)
        print(f'{gene}: T = {T:.2f} +- {np.sqrt(TV):.2f}, L = {L_rnd:.0f} +- {Lstd_rnd:.0f}')

        # Save to analyses
        analyses.loc[indices, 'T'] = T
        analyses.loc[indices, 'TV'] = TV
        analyses.loc[indices, 'L'] = L_est
        analyses.loc[indices, 'LV'] = LVest

    return analyses
