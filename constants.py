

import numpy as np

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


# %% Multithreading
cpu_count = 11

# %% Theory
slope_theory = 26   # pol/min
