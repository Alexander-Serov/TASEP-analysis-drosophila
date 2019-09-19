"""
Slope and max. polymerase number calculation procedure.
The characteristics are not calculated if less than 3 frames are present, i.e. nc duration >= 2*dt_new.

The function plots and outputs slope detection figures for each data set in the `output_slopes_folder`.

Also produces a histogram of recorded AP positions.
One AP per trace is used, and the AP is the mean observed AP value.
"""

import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

from constants import (AP_hist_folder, default_intensity_threshold, dt_new,
                       intensity_thresholds, output_slopes_folder,
                       slope_length_mins)
from Dataset import get_avg_on_regular_time_mesh
from set_figure_size import set_figure_size

# Plot paraemters
markersize = 2
lw = 1  # line width
colors = {'trace': '#68972F',
          'slope': '#ED1C24',
          'threshold': '#EDB120',
          'max': '#0072BD'}
alpha = 0.2     # transparency
height_factor = 0.5

ncs = range(11, 15)


def calculate_slopes(data, analyses_in, save_figures, pdf=False):
    """
    Main calculation procedure
    """
    analyses = analyses_in.copy()
    dataset_ids = set(data.dataset_id)

    for dataset_id in tqdm(dataset_ids, desc='Calculating slopes'):
        dataset_data = data[data.dataset_id == dataset_id]

        dataset_name = dataset_data.dataset.iloc[0]
        dataset_id = dataset_data.dataset_id.iloc[0]
        gene_id = dataset_data.gene_id.iloc[0]
        gene = dataset_data.gene.iloc[0]
        construct_id = dataset_data.construct_id.iloc[0]

        # Check if there are data points after filtering
        if not dataset_data.shape[0]:
            logging.warning(
                f"No data points available for slope fitting after positional fitering. Skipping dataset={dataset_id}")
            continue

        # %% Calculate slopes and max values
        set_figure_size(num=5, rows=1, page_width_frac=0.5,
                        clear=True, height_factor=height_factor)
        fig, ax = plt.subplots(1, 1, num=5)

        # Get an average trace on a regulare time mesh
        avg_data, std_data = get_avg_on_regular_time_mesh(dataset_data, dt=dt_new)

        # Plot the average trace
        x = avg_data.index
        y = avg_data.intensity
        ax.plot(x.values, y.values, '-o', markersize=markersize, lw=lw, c=colors['trace'])

        # Add standard deviation as a shaded region
        ylow = y - std_data.intensity
        yup = y + std_data.intensity
        plt.fill_between(x.values, ylow.values, yup.values,
                         facecolor=colors['trace'], alpha=alpha, edgecolor=None)

        # Detect the time limits of the plot
        t_span = avg_data[~pd.isna(avg_data.intensity)].index.values
        xlims = np.array(t_span[[0, -1]])
        # Extend slightly for better presentation
        xlims = xlims + (xlims[1] - xlims[0]) * 0.1 / 2 * np.array([-1, 1])
        plt.xlim(xlims)

        # %% Calculate slopes and max intensity
        for nc in ncs:
            nc_data = dataset_data[dataset_data.nc == nc]
            if not nc_data.empty:
                # Add columns to analyses
                cols = ['slope', 'slopeV', 'max', 'maxV']
                for col in cols:
                    if col not in analyses:
                        analyses[col] = np.nan

                # Fit slope on a fixed time interval after the start of the nc
                start_slope = nc_data.time.min()
                end_slope = start_slope + slope_length_mins

                # Calculate max and slope only if the nc includes at least 3 frames.
                # Otherwise, return nans
                nc_length = nc_data.time.max() - nc_data.time.min()
                if nc_length >= 2 * dt_new:
                    slope, slope_V, coefs = get_slope(nc_data, start=start_slope, end=end_slope)
                    max_intensity, max_intensity_std = get_max_intensity(nc_data)
                else:
                    slope, slope_V = [np.nan] * 2
                    coefs = [np.nan] * 2
                    max_intensity, max_intensity_std = [np.nan] * 2
                if slope < 0:
                    continue

                # Save to analyses
                analyses.loc[(dataset_id, nc), cols] = (
                    slope, slope_V, max_intensity, max_intensity_std**2)

                # Plot current nc slope and detected max inensity
                plot_slopes(dataset_data, coefs, start_slope, end_slope, max_intensity, nc, ax)

        # Save figure
        figname = "slopes_{id:02d}_{dataset_name}".format(id=dataset_id, dataset_name=dataset_name)
        figpath = os.path.join(output_slopes_folder, figname)

        if save_figures:
            if pdf:
                fig.savefig(figpath + '.pdf', pad_inches=0, bbox_inches='tight')

            # Save an enlarged png version for easier analysis
            factor = 3
            figsize = fig.get_size_inches()
            fig.set_figwidth(figsize[0] * factor)
            fig.set_figheight(figsize[1] * factor)
            fig.savefig(figpath + '.png', pad_inches=0, bbox_inches='tight')

        # %% Plot a histogram of recorded AP positions
        fig = plt.figure(num=10, clear=True)
        str_title = '%s, id=%i' % (dataset_name.replace('_', '\_'), dataset_id)

        MeanAPs = data.groupby(data['trace_id']).ap_mean.first()
        plt.hist(MeanAPs)

        # Adjust the plot
        plt.xlim([0, 1])
        plt.xlabel('AP')
        plt.title(str_title)

        # Save figure
        figpath = f"AP_hist_{dataset_id}.png"
        figpath = os.path.join(AP_hist_folder, figpath)
        if save_figures:
            fig.savefig(figpath)

    return analyses


def get_slope(nc_data, start, end):
    """
    Fit the first several minutes of a nuclear cycle with a straight line.
    Note that all points are fitted and not the average trace to avoid wrong weight assignment.
    A nc must have at least 4 points to be fitted.
    """
    cur_data = nc_data[(nc_data.time >= start)
                       & (nc_data.time <= end)]

    if len(cur_data) <= 3:
        slope = np.nan
        coefs = [np.nan, np.nan]
        slope_V = np.nan
    else:
        coefs, V = np.polyfit(cur_data.time, cur_data.intensity, 1, cov=True)
        slope = coefs[0]
        slope_V = V[0, 0]
    return slope, slope_V, coefs


def get_max_intensity(nc_data):
    """
    Get the maximum polymerase number for a cycle. Average within +- 1 time step of the regular mesh
    """

    if nc_data.empty:
        return np.nan

    # print(nc_data)
    avg_nc_data, _ = get_avg_on_regular_time_mesh(nc_data, dt_new)
    # print(avg_nc_data.intensity.dtype)
    max_time = avg_nc_data.intensity.idxmax()
    # print(max_time)

    # interval = max_time + np.array([-1, 1]) * dt_new
    int = nc_data[(nc_data.time >= max_time - dt_new) &
                  (nc_data.time <= max_time + dt_new)].intensity
    max_intensity = int.mean()
    max_intensity_std = np.std(int, ddof=1)
    return max_intensity, max_intensity_std


def plot_slopes(dataset_data, coefs, start, end, max_intensity, nc, ax):
    """
    Plot the initial slopes and measured max intensity of the nuclear cycle.
    Also plots the noise threshold.
    """
    dataset_name = dataset_data.dataset.iloc[0]
    dataset_id = dataset_data.dataset_id.iloc[0]

    # Plot the slope
    x_fit = np.asarray([start, end])
    y_fit = x_fit * coefs[0] + coefs[1]
    ax.plot(x_fit, y_fit, 'r-', lw=lw, c=colors['slope'])

    # Plot the noise threshold
    xlims = plt.xlim()
    intensity_threshold = intensity_thresholds.get(dataset_id, default_intensity_threshold)
    ax.plot(xlims, [intensity_threshold] * 2, 'g', lw=lw, c=colors['threshold'])

    # Print nc numbers over each detected nc
    _, ymax = plt.ylim()
    plt.text(start, ymax * 0.9, f'nc{nc}', fontsize=7)

    # Plot the max. polymerase number
    x_max = [dataset_data[dataset_data.nc == nc].time.min(),
             dataset_data[dataset_data.nc == nc].time.max()]
    y_max = np.asarray([1, 1]) * max_intensity
    plt.plot(x_max, y_max, 'm--', lw=lw, c=colors['max'])

    # Add plot title and adjust plot
    str_title = dataset_name.replace('_', '\_')
    plt.title(str_title)
    plt.xlabel('$t$, min')
    plt.ylabel('Fluo, a.u.')
    plt.ylim(ymin=-50)
