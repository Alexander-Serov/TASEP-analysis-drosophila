import numpy as np

from bayesian_linear_fit import bayesian_linear_fit
from constants import LV, L, k, kV
from least_squares_BCES import least_squares_BCES


def calculate_free_travel_time(analyses_in):
    """Estimate the free passage time by fitting the non-scaled current-density diagram across all ncs for each gene and construct individually
    """

    analyses = analyses_in.copy()
    genes = set(analyses.gene.values)
    constructs = set(analyses.construct.values)
    ncs = range(11, 15)
    columns = ['T', 'TV', 'L', 'LV']
    for column in columns:
        analyses[column] = np.nan

    # LV = 5e3 ** 2

    # Estimate the prior on the travel time
    T_prior = L / k
    TV_prior = T_prior**2 * (LV / L**2 + kV / k**2)

    def prior_func(m, sigmaV):
        return np.exp(-(m - T_prior)**2 / 2 / TV_prior) * (2 * np.pi * TV_prior)**(-1 / 2)

    # Group data by gene because must be same length for all constructs
    for gene_id, gene in enumerate(genes):
        # for construct_id, construct in enumerate(constructs):
        analyses_filtered = analyses[(analyses.gene == gene)]
        indices = analyses_filtered.index

        Ns, Ns_V, slopes, slopes_V = [], [], [], []
        # for nc in ncs:
        # Ns = np.concatenate([Ns, analyses_filtered['max' + str(nc)].values])
        Ns = analyses_filtered['max'].values
        NsV = analyses_filtered['maxV'].values
        slopes = analyses_filtered['slope'].values
        slopesV = analyses_filtered['slopeV'].values

        # Fit N/s passing through zero. Give result in minutes
        fit = bayesian_linear_fit(x=slopes, Vx=slopesV, y=Ns, Vy=NsV, c=0, prior=prior_func)
        T, TV = [fit[key] for key in ['m', 'mV']]

        # Estimate L from T and k
        L_est = k * T
        LVest = T**2 * kV + k**2 * TV
        # print(np.sqrt([kV / k**2, TV / T**2]))

        L_rnd, Lstd_rnd = np.round([L_est, np.sqrt(LVest)], decimals=-1)
        print(f'{gene}: T = {T:.2f} +- {np.sqrt(TV):.2f}, L = {L_rnd:.0f} +- {Lstd_rnd:.0f}')
        # print(fit)

        # Save
        analyses.loc[indices, 'T'] = T
        analyses.loc[indices, 'TV'] = TV
        analyses.loc[indices, 'L'] = L_est
        analyses.loc[indices, 'LV'] = LVest

    return analyses
