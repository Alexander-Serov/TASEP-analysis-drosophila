

import numpy as np
from numpy import sqrt

# %% Data import
data_folder = r"D:\Experimental_Data\Transcription. New data from Madhav (2016_07)"
matlab_csv_data_file = "all_data.csv"


# %% Slope detection
mins_per_frame = 37 / 60.
# expected_slope_length_frames = 5
detect_nc13_leftovers_interval_frames = 10
min_nc13_length_minutes = 7
max_nc13_length_minutes = 23
slope_length_frames = 7  # 7 frames corresponds to ~4 mins

# avg_slope_length_min = 3    # the length of the slope to look for with a sliding window
# avg_slope_length_frames = np.round(avg_slope_length_min / mins_per_frame).astype(np.int64)
# avg_slope_length_frames = 9  # the length of the slope to look for with a sliding window
# avg_slope_length_frames
output_slopes_folder = r".\output_slopes"
nc13_folder = r".\output_nc13"
ncs_locations_file = r".\ncs_locations.csv"
slopes_file = r".\slopes.csv"
max_reasonable_polymerase_number = 200
tex_slopes_data_file = r".\table_s_NSS_data.tex"


# %% Multithreading
cpu_count = 11

# %% TASEP Theory
k = 26 * 60  # polymerase elongation rate, pol/min
tau_rev = 5 / 60  # mean abortive initiation duration taken from Revyakin2006, min
l = 45  # polymerase footprint, bp
slope_simple_theory = k / (1 + sqrt(l))**2   # pol/min
slope_abortive_theory = (k * tau_rev - 1) / tau_rev / (k * tau_rev + l - 1)
slope_abortive_theory_error = 3   # pol/min
L_min = 5631  # minimum gene length (without 3' UTR), bp
N_simple_theory = L_min / sqrt(l) / (1 + sqrt(l))

# %% Gene expression patterns
hb_AP_interval = [0.25, 0.3]
kn_AP_interval = [0.6, 0.63]

# %% Bayes factors
n_pi = 4

# %% Misc
figures_folder = r".\figures"
