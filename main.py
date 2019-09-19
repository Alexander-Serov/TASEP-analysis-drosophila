"""
This file is the user interface for the analysis code in this folder. You can execute the code line-by-line to see the steps or execute it altogether to get the output figures produced in the `figures` folder.

Make sure  you deactivate the first `try... except...` block if you are launching the code from the terminal
"""

# Allow on-the-fly reloading of the modified project files for iPython-like environments. Disable if you are running in the console
try:
    has_run
except NameError:
    %matplotlib tk
    %load_ext autoreload
    %autoreload 2
    has_run = 1


import copy
# Package imports
import os

import matplotlib
import numpy as np
import pandas as pd
from tqdm import trange

from calculate_alpha import calculate_alpha
from calculate_free_travel_time import calculate_free_travel_time
from calculate_rho_and_J import calculate_rho_and_J
from calculate_slopes import calculate_slopes
from constants import (AP_hist_folder, data_folder, matlab_csv_data_file,
                       output_slopes_folder, s13_hb_bac)
# from Dataset import detect_timestep
from filter_by_AP import filter_by_AP
from identify_ncs import identify_ncs
from perform_welchs_test import perform_welchs_test
from plot_j_alpha_curve import plot_j_alpha_curve
from plot_normalized_current_density_diagram import \
    plot_normalized_current_density_diagram
from plot_parameter_evolution import plot_parameter_evolution
from reinit_folder import reinit_folder

# Do not show figures, just save to files
matplotlib.use('Agg')

# Constants
ncs = range(11, 15)

# %% Import
filepath = os.path.join(data_folder, matlab_csv_data_file)
data = pd.read_csv(filepath, sep=',', encoding='utf-8')
print('Loaded columns: ', data.columns.values)

# Calculate the mean AP position for each trace
data['ap_mean'] = data.ap.groupby(data['trace_id']).transform('mean')

# Detect the range of data sets and ncs present in the input file
datasets_len = data.dataset_id.max() + 1
datasets = set(data.dataset_id)
genes = set(data.gene)

# %% Perform AP filtering
filtered_data = filter_by_AP(data)

# %% Detect nuclear cycles
data_with_ncs, nc_limits = identify_ncs(filtered_data)
print('Time limits for each data set:\n', nc_limits)


# %% Initialize the data structure to store results
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


# %% Calculate initial slopes and maximum polymerase numbers
# You can modify the `save_figures` and `pdf` parameters if necessary

# Clean up the output folder
reinit_folder([AP_hist_folder, output_slopes_folder])
analyses = calculate_slopes(
    data_with_ncs, analyses, save_figures=True, pdf=False)


# %% Print the number of available data sets per gene, construct and nc
analyses['count'] = analyses.groupby(by=['gene', 'construct', 'nc']).Tstart.transform(
    'count')
print('Number of data sets per nc:\n', analyses.groupby(
    by=['gene', 'construct', 'nc']).Tstart.count())
print('\nNumber of data sets, all ncs mixed:\n', analyses.copy().reset_index().groupby(
    by=['gene', 'construct']).dataset_id.agg(lambda x: x.nunique()))


# %% Calculate the slope of the current-density diagram, the transit time and the effective length of the gene
analyses = calculate_free_travel_time(analyses)

# %% Get calibration coefficient I based on the nc13 measurements from (Zoller, Little, Gregor, 2018). Change `s13_hb_bac` in `constants.py` or `I_est` directly if you want to modify this behavior.
# Set recalibrate = False if you don't want to recalibrate the trace
recalibrate = False
if recalibrate:
    nc = 13
    avg_slope_13 = np.nanmean(
        analyses[(analyses.gene == 'hb') & (analyses.construct == 'bac')].loc[(slice(None), nc), :]['slope'].values)
    I_est = avg_slope_13 / s13_hb_bac    # fluo/pol
else:
    I_est = 25.28255294491828
print('Using the calibration coefficient I=', I_est)

# %% Use the calibration coeffeicient to estimate the polymerase flux, maximal number and alpha
analyses = calculate_rho_and_J(analyses, I_est)
analyses = calculate_alpha(analyses)


# %% Perform Welch's test for equal means
analyses = perform_welchs_test(analyses)


# %% Plot parameter evolution across ncs in different genes and constructs
for gene in genes:
    plot_parameter_evolution(analyses, gene=gene, pdf=True)


# %% Plot the current-density diagram and the dependencies of j and rho on alpha
analyses = plot_normalized_current_density_diagram(analyses, I_est, num=7, pdf=True)
plot_j_alpha_curve(analyses, pdf=True)

# %% Output all results into a .csv if necessary
# analyses.to_csv('analyses.csv')
