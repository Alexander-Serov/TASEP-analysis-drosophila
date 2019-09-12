
import matplotlib.pyplot as plt
import numpy as np

from constants import k, kV, l, lV, tau_rev, tau_revV

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
