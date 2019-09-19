
from constants import k, kV, l, lV
from theoretical_phase_parameters import J_over_k_MC
from theoretical_phase_parameters import rho_MC as rho_MC_func


def calculate_rho_and_J(analyses_in, I_est, IV_est=0):
    """
    The function uses the calibration coefficient to convert a.u. into polymerase numbers.
    It calculate:
    - the site occupation density rho,
    - the site density normalized to the maximal current regime r,
    - polymerase flux J/k in pol/min,
    - polymerase flux normalized to the MC regime

    One may provide the error on the calibration coefficient as `IV_est`
    """

    analyses = analyses_in.copy()
    T = analyses.T
    TV = analyses.TV

    rho_MC, rho_MCV = rho_MC_func()
    JoK_MC, JoK_MCV = J_over_k_MC()

    # rho
    max = analyses['max']
    maxV = analyses['maxV']
    analyses['rho'] = rho = max * l / I_est / k / T
    analyses['rhoV'] = rhoV = rho**2 * (maxV / max**2 + lV / l**2
                                        + IV_est / I_est**2 + kV / k**2 + TV / T**2)

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
