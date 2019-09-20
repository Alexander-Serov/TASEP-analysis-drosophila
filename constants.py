

import numpy as np
from numpy import sqrt

# %% Real data
# data_folder = r"D:\Experimental_Data\Transcription. Processed data from Ben"
# matlab_csv_data_file = "all_data.csv"

# Example data
data_folder = r".\example"
matlab_csv_data_file = "example_data.csv"


# avg_slope_length_min = 3    # the length of the slope to look for with a sliding window
# avg_slope_length_frames = np.round(avg_slope_length_min / mins_per_frame).astype(np.int64)
# avg_slope_length_frames = 9  # the length of the slope to look for with a sliding window
# avg_slope_length_frames
output_slopes_folder = r".\detected_slopes"
figures_folder = r'.\figure'
nc13_folder = r".\ncs_locations"
AP_hist_folder = r".\AP_hist"
tex_slopes_data_file = r".\table_s_NSS_data.tex"


# %% TASEP Theory
k = 26 * 60  # polymerase elongation rate, bp/min
kV = (2 * 60)**2  # var(k)
# k = 40 * 60  # polymerase elongation rate, bp/min, Fukaya2017 data
# kV = (5 * 60)**2  # var(k), Fukaya2017data

# tau_rev = 5 / 60  # mean abortive initiation duration taken from Revyakin2006, min
# tau_revV = (1 / 60)**2                    # var(tau_rev)
tau_rev = 3.4 / 60  # mean abortive initiation duration taken from Tantale2016, min
tau_rev * 60
tau_revV = (1 / 60)**2                    # var(tau_rev)

l = 50  # polymerase footprint, bp
lV = 8**2               # var(l)
slope_simple_theory = k / (1 + np.sqrt(l))**2   # pol/min
slope_abortive_theory = (k * tau_rev - 1) / tau_rev / (k * tau_rev + l - 1)
slope_abortive_theory_error = 3   # pol/min
max_pol_simple_theory = 109
max_pol_abortive_theory = max_pol_simple_theory * 0.30
L_min = 5631  # minimum gene length (without 3' UTR), bp
L = L_min
LV = 5e3**2             # var(L)
N_simple_theory = L_min / sqrt(l) / (1 + sqrt(l))

# %% Slope detection
mins_per_frame = 37 / 60.
min_nc13_length_minutes = 7
max_nc13_length_minutes = 23
# slope_length_frames = 6  # 7 frames corresponds to ~4 mins
slope_length_mins = L / k  # 4
slope_length_minsV = slope_length_mins**2 * (LV / L**2 + kV / k**2)  # 4
f'expected slope length = {slope_length_mins:.2f} +- {np.sqrt(slope_length_minsV):.2f}'


# %% Gene expression patterns
# hb_AP_interval = [0.25, 0.3]
# kn_AP_interval = [0.6, 0.64]

# %% Bayes factors
n_pi = 4

# %% Misc
figures_folder = r".\figures"

gene_labels = ['hb', 'kn', 'sn']
gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}
construct_labels = {0: 'b', 1: 'np', 2: 'ns'}

dt_new = 0.7  # min


# % Calibration
s13_hb_bac = 7  # pol/min calibration number from Ben's Cell 2018 article

# %% identify_ncs
# Noise thresholds
default_intensity_threshold = 40
intensity_thresholds = {0: 30,
                        1: 40,
                        2: 15,
                        3: 35,
                        5: 25,
                        6: 20,
                        7: 20,
                        8: 15,
                        9: 30,
                        10: 70,
                        11: 25,
                        14: 100,
                        33: 45,
                        34: 45,
                        39: 45,
                        40: 50,
                        41: 70,
                        44: 45,
                        45: 40,
                        47: 60,
                        50: 45,
                        51: 50,
                        52: 45,
                        53: 50,
                        56: 50,
                        58: 25,
                        59: 30,
                        60: 30,
                        64: 45,
                        66: 45, }

colors_additivity = {'bac': '#EDB120', 'no_pr': '#D95319',
                     'no_sh': '#68972F', 'sum': '#0072BD', 'abrt_limit': '#000000', 'brown': '#754C29'}
markers_additivity = ('x', 'o', '^', 'd', '*')

# Literature estimates
Darzacq2007 = {'label': 'Darzacq2007',
               'rho': 0.02282608695652174,
               'rhoV': 5.789224952741021e-05,
               'r': 0.026054183131503806,
               'rV': 9.071589255036988e-05,
               'j': 0.0008430185384405335,
               'jV': 9.292887120474494e-08,
               'alpha_over_k': 1.3223313433285069e-05,
               'alpha_over_kV': 1.864534134608124e-11,
               'tau': 30.667548190742508,
               'ktau': 2215.0108830550284}

Tantale2016 = {'label': 'Tantale2016',
               'rho': 0.19305019305019305,
               'rhoV': 0.02724142821376737,
               'JoK': 0.003861003861003861,
               'JoKV': 1.0514943104650967e-05,
               'r': 0.22035161317322577,
               'rV': 0.035342059231637614,
               'j': 0.2515140371572624,
               'jV': 0.04586314202098771,
               'alpha_over_k': 0.004761904761904762,
               'alpha_over_kV': 2.5787241415756705e-05,
               'tau': 3.4054054054054053,
               'tauV': 16.92377297812088}


1 - 1 / (1 + np.sqrt(45))
L / k
1 / (1 + np.sqrt(12))
1 / k / tau_rev

k * np.array([3.65, 4.23, 3.99])

np.array([2.7, 0.3]) * 1000 / 60
