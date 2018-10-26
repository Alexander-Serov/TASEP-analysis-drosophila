"""
TODO:
"""

try:
    has_run
except NameError:
    %matplotlib tk
    %load_ext autoreload
    %autoreload 2
    has_run = 1
else:
    print("Graphic interface NOT re-initialized")


import os
# from find_best_fit import find_best_fit
from functools import partial
from multiprocessing import Pool

# from identify_shift_and_slope import identify_shift_and_slope
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import trange

from calculate_bayes_factor import calculate_bayes_factor
from calculate_slopes import calculate_slopes
from constants import (AP_hist_folder, cpu_count, data_folder,
                       detect_nc13_leftovers_interval_frames, gene_labels,
                       matlab_csv_data_file, max_nc13_length_minutes,
                       min_nc13_length_minutes, mins_per_frame, n_pi,
                       nc13_folder, ncs_locations_file, output_slopes_folder,
                       slope_length_frames, slopes_file)
from detect_y_regions_in_snail import detect_y_regions_in_snail
from identify_ncs import identify_ncs
from plot_boxplot_hb import plot_boxplot_hb
from plot_max_histograms import plot_max_histograms
from plot_max_num_boxplot import plot_max_num_boxplot
from plot_one_embryo import plot_one_embryo, plot_xy_in_all_embryos
from plot_slope_histograms import plot_slope_histograms
from plot_slopes_boxplot import plot_slopes_boxplot
from reinit_folder import reinit_folder
from save_series_plot import save_series_plot
from slopes_to_tex import slopes_to_tex

# # %% Initialize


# %% Import
filepath = os.path.join(data_folder, matlab_csv_data_file)
data = pd.read_csv(filepath, sep=';')

# Calculate mean xPos and yPos for each particle
data['xPosMean'] = data['xPos'].groupby(data['trace_id']).transform('mean')
data['yPosMean'] = data['yPos'].groupby(data['trace_id']).transform('mean')

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
    # reinit_folder(output_slopes_folder)
    reinit_folder(nc13_folder)
    ncs_locations = identify_ncs(data)
    # Save
    ncs_locations.to_csv(filepath, sep=';')
ncs_locations

# %%
# data[data.dataset_id == 8]

# %% Calculate the slopes by group
filepath = os.path.join(data_folder, slopes_file)
bl_load = True
# bl_load = False
if bl_load and os.path.exists(filepath):
    slopes = pd.read_csv(filepath, sep=';')
else:
    # Initialize
    slopes = pd.DataFrame(columns=['dataset_id', 'gene_id', 'construct_id', 'AP',
                                   'slope_nc13', 'slope_nc14', 'max_nc13', 'max_nc14', 'slope_nc13_count', 'slope_nc14_count'], index=range(datasets_len), dtype=np.float64)

    slope_start_pols = []
    for dataset_id in trange(datasets_len):
        slopes.iloc[dataset_id], slope_start_pol = calculate_slopes(
            data[data.dataset_id == dataset_id], ncs_locations)
        slope_start_pols += slope_start_pol

    # Save
    slopes.to_csv(filepath, sep=';')

slopes.dtypes
slopes
slopes.groupby(by=['gene_id', 'construct_id']).count()
# %%
# srt_lst = sorted([pol for pol in slope_start_pols if not np.isnan(pol)])

# reinit_folder(r'.\xy_plots')
# plot_xy_in_all_embryos(data)

# %%
# calculate_slopes(data[data.dataset_id == 28], ncs_locations)

# %% Detect y regions in the snail data
# detect_y_regions_in_snail(data)


# # %% Plot the histograms of the slopes distribution
# plot_slope_histograms(slopes)
# plot_max_histograms(slopes)

# %% Boxplots
# plot_slopes_boxplot(slopes)
# plot_max_num_boxplot(slopes)
plot_boxplot_hb(slopes)


# %% Calculate Bayes factors for slopes

# All-slopes parameters for the prior
all_slopes = np.concatenate([slopes.slope_nc13.dropna(), slopes.slope_nc14.dropna()])
mu_pi = all_slopes.mean()
V_pi = all_slopes.var()


def combination(gene_id, construct_id, nc):
    return 6 * gene_id + 2 * construct_id + (nc - 13)


construct_labels = ['bac', 'no_pr', 'no_sh']
labels = []
for gene_id in range(3):
    for construct_id in range(3):
        for nc in [13, 14]:
            label = '%s_%s_nc_%i' % (gene_labels[gene_id], construct_labels[construct_id], nc)
            labels.append(label)
labels


nc_labels = ['slope_nc13', 'slope_nc14']

bayes_factors = np.ones([18, 18], dtype=np.float64) * np.nan

for gene_id1 in range(3):
    for construct_id1 in range(3):
        for gene_id2 in range(3):
            for construct_id2 in range(3):
                for nc1 in [13, 14]:
                    nc_label1 = nc_labels[nc1 - 13]
                    for nc2 in [13, 14]:
                        nc_label2 = nc_labels[nc2 - 13]
                        # 13 - 13
                        dataset1 = slopes[(slopes.gene_id == gene_id1) & (
                            slopes.construct_id == construct_id1)][nc_label1].values
                        dataset2 = slopes[(slopes.gene_id == gene_id2) & (
                            slopes.construct_id == construct_id2)][nc_label2].values
                        bayes_factors[combination(gene_id1, construct_id1, nc1), combination(
                            gene_id2, construct_id2, nc2)], _, _ = calculate_bayes_factor(dataset1, dataset2, n_pi, mu_pi, V_pi)

# Set diagonal to nan
for i in range(18):
    bayes_factors[i, i] = np.nan

# Threshold
B_log_lim = 1
B_trsh = bayes_factors * 0
B_trsh[bayes_factors >= B_log_lim] = 1
B_trsh[bayes_factors <= -B_log_lim] = -1

# Mask
for i in range(18):
    for j in range(18):
        if j >= i:
            B_trsh[i, j] = np.nan

bayes_factors
plt.rc('text', usetex=False)
fig = plt.figure(num=4)
fig.clf()
fig, ax = plt.subplots(num=4)
# ax.cla()
# fig.clf()
im = ax.imshow(B_trsh[0:6, 0:6], cmap='inferno')

# ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
# Show all ticks
ax.set_xticks(np.arange(len(labels) / 3))
ax.set_yticks(np.arange(len(labels) / 3))
# set labels
ax.set_xticklabels(labels[0:6])
ax.set_yticklabels(labels[0:6])
# Rotate Labels
# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

ax.figure.colorbar(im, ax=ax)
fig.tight_layout()
plt.show()


#

# %% Output slopes table to tex file
slopes_to_tex(slopes)


# %%

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
