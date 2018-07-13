

import numpy as np

# %% Data import
matlab_csv_filepath = r"D:\Experimental_Data\Transcription. New data from Madhav (2016_07)\all_data.csv"


# %% Slope detection
mins_per_frame = 37 / 60.
expected_slope_length_frames = 5
detect_nc13_leftovers_interval_frames = 10
# avg_slope_length_min = 3    # the length of the slope to look for with a sliding window
# avg_slope_length_frames = np.round(avg_slope_length_min / mins_per_frame).astype(np.int64)
# avg_slope_length_frames = 9  # the length of the slope to look for with a sliding window
# avg_slope_length_frames
output_slopes_folder = r".\output_slopes"
nc13_folder = r".\output_nc13"
