

try:
    has_run
except NameError:
    %matplotlib tk
    %load_ext autoreload
    %autoreload 2
    has_run = 1
else:
    print("Graphic interface NOT re-initialized")


from calculate_slopes import calculate_slopes
# from find_best_fit import find_best_fit
from functools import partial
from identify_ncs import identify_ncs
# from identify_shift_and_slope import identify_shift_and_slope
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
from plot_slope_histograms import plot_slope_histograms
from reinit_folder import reinit_folder
from save_series_plot import save_series_plot
from tqdm import trange

from constants import data_folder, matlab_csv_data_file, mins_per_frame, output_slopes_folder, detect_nc13_leftovers_interval_frames, nc13_folder, min_nc13_length_minutes, max_nc13_length_minutes, slope_length_frames, cpu_count, ncs_locations_file, slopes_file


# # %% Initialize
# reinit_folder(output_slopes_folder)
# reinit_folder(nc13_folder)


# %% Import
filepath = os.path.join(data_folder, matlab_csv_data_file)
data = pd.read_csv(filepath, sep=';')
print(data.columns.values)

# %% Analyze
traces_len = data.trace_id.max() + 1
datasets_len = data.dataset_id.max() + 1
intersects = np.ones([traces_len, 1]) * np.nan
slopes = np.ones([traces_len, 1]) * np.nan

# %% Identify the locations of nc13 and nc 14 in the data
bl_load = True
# bl_load = False
filepath = os.path.join(data_folder, ncs_locations_file)
if bl_load & os.path.exists(filepath):
    ncs_locations = pd.read_csv(filepath, sep=';')
else:
    # Initialize folders
    reinit_folder(output_slopes_folder)
    reinit_folder(nc13_folder)
    ncs_locations = identify_ncs(data)
    # Save
    ncs_locations.to_csv(filepath, sep=';')
# ncs_locations


# %% Calculate the slopes in each trace
bl_load = True
# bl_load = False
filepath = os.path.join(data_folder, slopes_file)
if bl_load and os.path.exists(filepath):
    slopes = pd.read_csv(filepath, sep=';')
else:
    # Initialize
    slopes = pd.DataFrame(columns=['dataset_id', 'gene_id', 'construct_id',
                                   'slope_nc13', 'slope_nc14'], index=range(traces_len))

    for trace_id in trange(0 + 1 * traces_len):
        slopes.iloc[trace_id] = calculate_slopes(
            data[data.trace_id == trace_id], ncs_locations, trace_id)

    # Save
    slopes.to_csv(filepath, sep=';')

# slopes


# %% Plot the histograms of the slopes distribution
plot_slope_histograms(slopes)

#

#


# cur_slopes

# slopes.slope_nc13.hist()

#

#
# pool = Pool(cpu_count)
#
#
# # def short_func(trace_id):
#     return calculate_slopes(data[data.trace_id == trace_id], ncs_locations, trace_id)
#
#
# slopes_list = pool.map(short_func, trange(10 + 0 * traces_len))
#
# # for trace_id in trange(1000 + 0 * traces_len):
# # slopes.iloc[trace_id] = short_func(trace_id)
# # # Save
# # slopes.iloc[trace_id] = pd.Series({'trace_id': trace_id, 'dataset_id': dataset_id, 'gene_id':
# #                                    gene_id, 'construct_id': construct_id, 'slope_nc13': slope_nc13, 'slope_nc14': slope_nc14})
#
#
# # pool.close()
# slopes_list

#
