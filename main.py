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


import itertools
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

from calculate_alpha import calculate_alpha
from calculate_bayes_factor import calculate_bayes_factor
from calculate_free_travel_time import calculate_free_travel_time
from calculate_rho_and_J import calculate_rho_and_J
from calculate_slopes import calculate_slopes
from constants import (AP_hist_folder, L, cpu_count, data_folder,
                       detect_nc13_leftovers_interval_frames, dt_new,
                       gene_labels, k, kV, l, lV, matlab_csv_data_file,
                       max_nc13_length_minutes, max_pol_abortive_theory,
                       max_pol_simple_theory, min_nc13_length_minutes,
                       mins_per_frame, n_pi, nc13_folder, ncs_locations_file,
                       output_slopes_folder, s13_hb_bac, slope_abortive_theory,
                       slope_simple_theory, slopes_file, tau_rev)
from Dataset import detect_timestep, get_regular_time_mesh
from detect_y_regions_in_snail import detect_y_regions_in_snail
from filter_by_AP import filter_by_AP
from identify_ncs import identify_ncs
from perform_additivity_test import perform_additivity_test
# % Additivity plots
from plot_additivity_plots import plot_additivity_plots
from plot_alpha import plot_alpha
from plot_boxplot_hb import plot_boxplot_hb
from plot_current_density_diagram import plot_current_density_diagram
from plot_j_alpha_curve import plot_j_alpha_curve
from plot_max_histograms import plot_max_histograms
from plot_max_num_boxplot import plot_max_num_boxplot
from plot_normalized_current_density_diagram import \
    plot_normalized_current_density_diagram
from plot_one_embryo import plot_one_embryo, plot_xy_in_all_embryos
from plot_slope_histograms import plot_slope_histograms
from plot_slopes_boxplot import plot_slopes_boxplot
# Add theory
from plot_theory_current_density import plot_theoretical_curve
from reinit_folder import reinit_folder
from save_series_plot import save_series_plot
# % Current vs density plot
from set_figure_size import set_figure_size
from slopes_to_tex import slopes_to_tex
# Abortive prediction
from theoretical_phase_parameters import (J_over_k_LD, alpha_over_k_abortive,
                                          alpha_over_k_MC, rho_LD)

# Constants
ncs = range(11, 15)


# %% Import
filepath = os.path.join(data_folder, matlab_csv_data_file)
data = pd.read_csv(filepath, sep=';')

dt, _, _ = detect_timestep(data)


# Calculate mean xPos and yPos for each particle
# data['xPosMean'] = data['xPos'].groupby(data['trace_id']).transform('mean')
# data['yPosMean'] = data['yPos'].groupby(data['trace_id']).transform('mean')
data['ap_mean'] = data.ap.groupby(data['trace_id']).transform('mean')
# data['ap_median'] = data.ap.groupby(data['trace_id']).transform('median')

print(data.columns.values)

# # Create a regular time mesh
# time_mesh_regular = get_regular_time_mesh(data, dt_new)


# %% Analyze
traces_len = data.trace_id.max() + 1
datasets_len = data.dataset_id.max() + 1
intersects = np.full([traces_len, 1], np.nan)
slopes = np.full([traces_len, 1], np.nan)
datasets = set(data.dataset_id)

# %% #TODO AP filtering here
data = filter_by_AP(data)

# data[data.dataset_id == 13]


# # %% Identify the locations of nc13 and nc 14 in the data
data_with_ncs, nc_limits = identify_ncs(data)
nc_limits.loc[13, :]


# Calculate the slopes for each dataset (in arbitrary units)
# Prepare the analyses object
genes = set(data.gene.values)
constructs = set(data.construct.values)
index = pd.MultiIndex.from_product((datasets, ncs), names=['dataset_id', 'nc'])

analyses = pd.DataFrame(columns=['dataset_name', 'gene', 'gene_id',
                                 'construct'], index=index)
# Copy dataset names, genes and constructs
for id in datasets:
    analyses.loc[id, 'gene'] = data[data.dataset_id == id].head(1).gene.values[0]
    analyses.loc[id, 'gene_id'] = data[data.dataset_id == id].head(1).gene_id.values[0]
    analyses.loc[id, 'construct'] = data[data.dataset_id == id].head(1).construct.values[0]
    analyses.loc[id, 'dataset_name'] = data[data.dataset_id == id].head(1).dataset.values[0]

analyses = pd.concat([analyses, nc_limits], axis='columns')


filepath = os.path.join(data_folder, slopes_file)
# bl_load = True
bl_load = False
save_figures = 1
if bl_load and os.path.exists(filepath):
    slopes = pd.read_csv(filepath, sep=';')
else:
    # Initialize
    # slopes = pd.DataFrame(columns=['dataset_id', 'gene_id', 'construct_id', 'AP',
    #                                'slope_nc13', 'slope_nc14', 'max_nc13', 'max_nc14', 'slope_nc13_count', 'slope_nc14_count'], index=range(datasets_len), dtype=np.float64)
    reinit_folder([AP_hist_folder, output_slopes_folder])

    slope_start_pols = []
    for dataset_id in trange(datasets_len):
        # dataset_id = 2
        # print(dataset_id)
        analyses = calculate_slopes(
            data_with_ncs[data_with_ncs.dataset_id == dataset_id], analyses, save_figures=save_figures, pdf=0)
    # break

# print(analyses)

# %% Calculate the number of data points
analyses['count'] = analyses.groupby(by=['gene', 'construct', 'nc']).Tstart.transform(
    'count')  # .rename('count').reset_index()
print(analyses.groupby(by=['gene', 'construct', 'nc']).Tstart.count())
print(analyses.copy().reset_index().groupby(
    by=['gene', 'construct']).dataset_id.agg(lambda x: x.nunique()))


# %%
analyses = calculate_free_travel_time(analyses)
# analyses.loc[:, ['construct', 'L', 'LV']]

# %% Get calibration I based on Ben's nc13 measurement
nc = 13
avg_slope_13 = np.nanmean(
    analyses[(analyses.gene == 'hb') & (analyses.construct == 'bac')].loc[(slice(None), nc), :]['slope'].values)
I_est = avg_slope_13 / s13_hb_bac    # fluo/pol
# print(I_est)

analyses = calculate_rho_and_J(analyses, I_est)
analyses[analyses.gene == 'kn'][['construct', 'rho', 'rhoV', 'JoK', 'JoKV']]
analyses = calculate_alpha(analyses, I_est)
# analyses =
# analyses[analyses.gene == 'kn'][['construct', 'rho', 'rhoV',
#                                  'JoK', 'JoKV', 'alpha_over_k_J', 'alpha_over_k_JV', 'alpha_over_k_rho', 'alpha_over_k_rhoV',  'alpha_over_k_comb', 'alpha_over_k_combV']]

# %%
# analyses.columns
# Index(['dataset_name', 'gene', 'construct', 'Tstart', 'Tend', 'slope',
#        'slopeV', 'max', 'maxV', 'count', 'T', 'TV', 'L', 'LV', 'rho', 'rhoV',
#        'JoK', 'JoKV', 'alpha_over_k_J', 'alpha_over_k_rho', 'kappa',
#        'alpha_origin', 'alpha_over_k_comb', 'alpha_over_k_combV', 'tau',
#        'tauV', 'rho_p_value', 'JoK_p_value', 'alpha_over_k_comb_p_value',
#        'tau_p_value'],
#       dtype='object')
# analyses.loc[:, ('construct', 'kappa', 'alpha_origin', 'alpha_over_k_comb', 'alpha_over_k_combV')]
analyses = perform_additivity_test(analyses)
analyses[['rho', 'j']].max()
# analyses.columns
# analyses.loc[(slice(None), 13), ['gene', 'abortive_theory_p_value']].groupby('gene').first()


# print(Ts)


# %% Make additivity plots
for gene in 'hb kn sn'.split():
    plot_additivity_plots(analyses, gene=gene, pdf=1)

analyses.groupby(by=['gene', 'nc']).first()
# construct_avg = analyses_norm.groupby(by=['gene', 'construct']).mean()
# construct_std = analyses_norm.groupby(by=['gene', 'construct']).std()


# %% Calculate alpha. Automatically saved to analyses. Plot then the evolution of alpha acros ncs and genes

# print(analyses.index)
# print(analyses)
analyses.to_csv('analyses.csv', sep=';')
# plot_alpha(analyses)


# %%
print('Number of identified slopes per group:')
print(analyses.copy().reset_index().groupby(
    by=['gene', 'construct', 'nc']).JoK.agg(lambda x: x.nunique()).astype(int))
print(analyses.columns.values)

# np.array([0.11, np.sqrt(0.000663)]) * 60

# %% Numerical estimates for the article
r1 = analyses.copy().reset_index().groupby(
    by=['gene', 'construct', 'nc']).first().loc[:, ['JoK', 'JoKV', 'rho', 'rhoV']]  # .agg(lambda x: x.nunique()))
r1['JoKstd'] = r1.agg(lambda x: np.sqrt(x.JoKV), axis=1)
r1['rhostd'] = r1.agg(lambda x: np.sqrt(x.rhoV), axis=1)
print(r1)
print(r1 / J_over_k_LD(alpha_over_k_MC) * 100)


# print(f'J/k = {r1a.JoK} +- {np.sqrt(r1a.JoKV)}')

# %%
# analyses = plot_current_density_diagram(analyses, I_est, num=7, pdf=1)
analyses = plot_normalized_current_density_diagram(analyses, I_est, num=7, pdf=1)
# %%
plot_j_alpha_curve(analyses, pdf=1)
# analyses.columns

# %% Back-of-the-envelope estimates for the critical alpha

1 / (1 + np.sqrt(50))
8 / (2 * np.sqrt(50) * (1 + np.sqrt(50))**2)

part_k = 1 / (1 + np.sqrt(l))
part_l = -k / 2 / np.sqrt(l) / (1 + np.sqrt(l))**2
a_crit = k / (1 + np.sqrt(l))   # in min^-1
a_critV = np.sqrt(kV * part_k**2 + lV * part_l**2)

print(f'alpha^* = {a_crit} +- {np.sqrt(a_critV)}')

# dalpha_crit = np.sqrt(
