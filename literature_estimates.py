import numpy as np
from scipy.special import erf, gammaln

from BLUE_estimator import BLUE_estimator
from constants import k as k_Bothma

# % Common constants
l = 50
JoK_max = (1 + np.sqrt(l))**(-2)

nmRNA = np.array([200, 400])

300 * 0.3 * 50 / 2300
t1 = np.array([200, 400]) * 0.7
t1 * 50 / 2300

np.array([1, 2]) * 50 * 0.7 / 2300

J1 = nmRNA * 0.39 * 0.00159 / 200
kElong = 2300 / 32
print(J1 / kElong)

50 / (4.1 * 4000 / 60)
1 * 60 / 4.1 / 4000


def alpha_over_k_rho_func(rho, rhoV, l, lV):
    """Inverse function for alpha/k from J/k"""
    aoK = rho / (l - rho * (l - 1))
    partial_rho = aoK**2 * l / rho**2
    partial_l = - aoK**2 * (1 - rho) / rho

    aoKV = partial_rho**2 * rhoV + partial_l**2 * lV

    return [aoK, aoKV]


def alpha_over_k_J_func(JoK, JoKV, l, lV):
    """Inverse function for alpha/k from J/k"""
    aoK = (1
           - JoK * (l - 1)
           - np.sqrt(1 + JoK**2 * (l - 1) **
                     2 - 2 * JoK * (l + 1))
           ) / 2
    partial_l = JoK / 2 * (-1 + (1 + JoK - JoK * l) / np.sqrt(-4 * JoK + (-1 + JoK * (l - 1))**2))
    partial_JoK = (1 - l + (1 - JoK * (l - 1)**2 + l) /
                   np.sqrt(-4 * JoK + (-1 + JoK * (l - 1))**2)) / 2
    aoKV = partial_l**2 * lV + partial_JoK**2 * JoKV
    return [aoK, aoKV]


# %% Tantale2016 estimates
l = 50
lV = 8**2
dt = 4.2  # s
dtV = 2.6**2
k = 3.7 * 1000 / 60  # bp/s
kV = (2.1 * 1000 / 60)**2
Tantale2016 = {'label': 'Tantale2016'}

Tantale2016['rho'] = rho = l / dt / k
Tantale2016['rhoV'] = rhoV = rho**2 * (lV / l**2 + dtV / dt**2 + kV / k**2)
print('rho', [rho, np.sqrt(rhoV)])

Tantale2016['JoK'] = JoK = 1 / dt / k
Tantale2016['JoKV'] = JoKV = JoK**2 * (dtV / dt**2 + kV / k**2)
print('JoK', [JoK, np.sqrt(JoKV)])
print('J / Jmax', JoK / JoK_max)

# Normalized parameters:
Tantale2016['r'] = r = (l + np.sqrt(l)) / k / dt
partial_l = (1 + 1 / 2 / np.sqrt(l)) / k / dt
partial_k = - r / k
partial_dt = - r / dt
Tantale2016['rV'] = rV = partial_l**2 * lV + partial_k**2 * kV + partial_dt**2 * dtV
print('r', [r, np.sqrt(rV)])

Tantale2016['j'] = j = (1 + np.sqrt(l))**2 / k / dt
partial_l = (1 + np.sqrt(l)) / k / dt / np.sqrt(l)
partial_k = - j / k
partial_dt = - j / dt
Tantale2016['jV'] = jV = partial_l**2 * lV + partial_k**2 * kV + partial_dt**2 * dtV
print('j', [j, np.sqrt(jV)])

alpha_over_k, alpha_over_kV = alpha_over_k_rho_func(rho, rhoV, l, lV)
Tantale2016['alpha_over_k'] = alpha_over_k
Tantale2016['alpha_over_kV'] = alpha_over_kV
print('alpha_over_k', [alpha_over_k, np.sqrt(alpha_over_kV)])

Tantale2016['tau'] = tau = 1 / k / alpha_over_k
Tantale2016['tauV'] = tauV = tau**2 * (kV / k**2 + alpha_over_kV / alpha_over_k**2)
print('tau', [tau, np.sqrt(tauV)])

# ktau estimate
ktau = 1 / alpha_over_k
ktauV = ktau**2 * (alpha_over_kV / alpha_over_k**2)
print('ktau', [ktau, np.sqrt(ktauV)])
print(Tantale2016)
Tantale2016

# %% Darzacz2007 estimates
l = 50
lV = 8**2
nGenes = 200
L = 2300
nMRNA = 300
nMRNAV = 100 ** 2
RMS2 = 0.7
T = 32  # s
RInit = 0.39
kEscape = 0.00159   # s^{-1}
Darzacq2007 = {'label': 'Darzacq2007'}


Darzacq2007['rho'] = rho = nMRNA * RMS2 * l / L / nGenes
Darzacq2007['rhoV'] = rhoV = rho**2 * nMRNAV / nMRNA**2
print('rho', [rho, np.sqrt(rhoV)])

k = L / T
k * 60 / 1e3
J = nMRNA * RInit * kEscape / nGenes
JoK = J / k
JoKV = JoK**2 * nMRNAV / nMRNA**2
print('JoK', [JoK, np.sqrt(JoKV)])
print('J / Jmax', JoK / JoK_max)

# Normalized parameters:
r = nMRNA * RMS2 * (l + np.sqrt(l)) / L / nGenes
partial_l = r / (l + np.sqrt(l)) * (1 + 1 / 2 / np.sqrt(l))
partial_nMRNA = r / nMRNA
rV = partial_l**2 * lV + partial_nMRNA**2 * nMRNAV
Darzacq2007['r'] = r
Darzacq2007['rV'] = rV
print('r', [r, np.sqrt(rV)])

j = nMRNA * RInit * kEscape / nGenes / k * (1 + np.sqrt(l))**2
partial_l = j / (l + np.sqrt(l))
partial_nMRNA = j / nMRNA
jV = partial_l**2 * lV + partial_nMRNA**2 * nMRNAV
Darzacq2007['j'] = j
Darzacq2007['jV'] = jV
print('j', [j, np.sqrt(jV)])


alpha_over_k_rho, alpha_over_k_rhoV = alpha_over_k_rho_func(rho, rhoV, l, lV)
print('alpha_over_k_rho', [alpha_over_k_rho, np.sqrt(alpha_over_k_rhoV)])

alpha_over_k_J, alpha_over_k_JV = alpha_over_k_J_func(JoK, JoKV, l, lV)
print('alpha_over_k_J', [alpha_over_k_J, np.sqrt(alpha_over_k_JV)])

# Get a blue estimator for alpha_over_k
alpha_over_k, alpha_over_kV, kappa = BLUE_estimator(
    alpha_over_k_rho, alpha_over_k_rhoV, alpha_over_k_J, alpha_over_k_JV)
Darzacq2007['alpha_over_k'] = alpha_over_k
Darzacq2007['alpha_over_kV'] = alpha_over_kV
print('alpha_over_k', [alpha_over_k, np.sqrt(alpha_over_kV)])

tau_rho = 1 / k / alpha_over_k_rho
tau_rhoV = tau_rho**2 * (kV / k**2 + alpha_over_k_rhoV / alpha_over_k_rho**2)
print('tau_rho', [tau_rho, np.sqrt(tau_rhoV)])

tau_J = 1 / k / alpha_over_k_J
tau_JV = tau_J**2 * (kV / k**2 + alpha_over_k_JV / alpha_over_k_J**2)
print('tau_J', [tau_J, np.sqrt(tau_JV)])

# Get a blue estimator for tau
tau, tauV, kappa = BLUE_estimator(tau_rho, tau_rhoV, tau_J, tau_JV)
Darzacq2007['tau'] = tau
print('tau', [tau, np.sqrt(tauV)])

# Get a blue estimator for k*tau
ktau_rho = 1 / alpha_over_k_rho
ktau_rhoV = ktau_rho**2 * (alpha_over_k_rhoV / alpha_over_k_rho**2)
print('ktau_rho', [ktau_rho, np.sqrt(ktau_rhoV)])

ktau_J = 1 / alpha_over_k_J
ktau_JV = ktau_J**2 * (alpha_over_k_JV / alpha_over_k_J**2)
print('ktau_J', [ktau_J, np.sqrt(ktau_JV)])

ktau, ktauV, kappa = BLUE_estimator(ktau_rho, ktau_rhoV, ktau_J, ktau_JV)
Darzacq2007['ktau'] = ktau
print('ktau', [ktau, np.sqrt(ktauV)])
print(Darzacq2007)
Darzacq2007


# %% non-related
n = 3.5e5
ln_ptie = - gammaln(n / 2 + 1) * 2 + gammaln(n + 1) - (n - 1) * np.log(2)
np.exp(ln_ptie)
lg_ptie = ln_ptie / np.log(10)
p_tieall = 1 - (1 - np.exp(ln_ptie)) ** 435
20597 / 435 * 2

# %%
# k is in bp/min
1 / (0.03 * k_Bothma) * 60
