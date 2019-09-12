import numpy as np

from constants import k, kV, l, lV
from least_squares_BCES import least_squares_BCES
from theoretical_phase_parameters import J_over_k_MC
from theoretical_phase_parameters import rho_MC as rho_MC_func


def calculate_rho_and_J(analyses_in, I_est, IV_est=0):
    """
    """

    analyses = analyses_in.copy()
    genes = set(analyses.gene)
    constructs = set(analyses.construct)
    ncs = set(analyses.index.get_level_values('nc'))
    # columns = ['rho', 'rhoV', 'JoK', 'JoKV']
    # for column in columns:
    #     analyses[column] = np.nan
    T = analyses['T']
    TV = analyses.TV

    rho_MC, rho_MCV = rho_MC_func()
    JoK_MC, JoK_MCV = J_over_k_MC()

    # cols = ['rho_comb', 'rho_combV', 'JoK_comb', 'JoK_combV', 'rho_JoK_comb_cov']
    # for col in cols:
    #     analyses[col] = np.nan

    # for nc in ncs:
    # rho
    max = analyses['max']
    maxV = analyses['maxV']
    analyses['rho'] = rho = max * l / I_est / k / T
    analyses['rhoV'] = rhoV = rho**2 * (maxV / max**2 + lV / l**2 + IV_est /
                                        I_est**2 + kV / k**2 + TV / T**2)

    analyses['r'] = r = rho / rho_MC
    analyses['rV'] = r**2 * (rhoV / rho**2 + rho_MCV / rho_MC**2)

    # J/k
    slope = analyses['slope']
    slopeV = analyses['slopeV']
    analyses['JoK'] = JoK = slope / k / I_est
    analyses['JoKV'] = JoKV = JoK**2 * (slopeV / slope**2 + kV / k**2 + IV_est / I_est**2)

    analyses['j'] = j = JoK / JoK_MC
    analyses['jV'] = j**2 * (JoKV / JoK**2 + JoK_MCV / JoK_MC**2)

    return analyses
