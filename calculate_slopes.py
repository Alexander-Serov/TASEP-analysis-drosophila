

import logging
import os
import sys
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constants import (AP_hist_folder, L, default_intensity_threshold, dt_new,
                       intensity_thresholds, k, mins_per_frame,
                       output_slopes_folder, slope_length_mins)
from Dataset import get_avg_on_regular_time_mesh
from set_figure_size import set_figure_size

markersize = 2
lw = 1
colors = {'trace': '#68972F',
          'slope': '#ED1C24',
          'threshold': '#EDB120',
          'max': '#0072BD'}
alpha = 0.2


def calculate_slopes(data, analyses_in, save_figures, pdf=False):
    """
    A new slope calculation procedure for Ben's data.
    Calculate nc13 and nc14 slope for the given dataset.
    Calculate max polymerase numbers in these nc as well.

    The characteristics are not calculated if less than 3 frames are present, i.e. nc duration >= 2*dt_new
    """
    analyses = analyses_in.copy()

    dataset_name = data.dataset.iloc[0]
    dataset_id = data.dataset_id.iloc[0]
    gene_id = data.gene_id.iloc[0]
    gene = data.gene.iloc[0]
    construct_id = data.construct_id.iloc[0]

    # AP filtering (one AP per trace)
    MeanAPs = data.groupby(data['trace_id']).ap_mean.first()
    # max_AP = MeanAPs.max()
    # median_AP = MeanAPs.median()
    height_factor = 0.5

    gene_data = data
    # if gene == 'hb':
    #     # do not include the transition region of width 0.05 at the right side
    #     AP_trans_region_width = 0.05
    #     gene_data = gene_data[(gene_data.ap_mean <= max_AP - AP_trans_region_width)
    #                           & (gene_data.ap_mean >= 0.05)]
    #     # print(max_AP)
    # elif gene == 'kn':
    #     AP_central_region_width = 0.04
    #     gene_data = gene_data[(gene_data.ap_mean >= median_AP - AP_central_region_width / 2) &
    #                           (gene_data.ap_mean <= median_AP + AP_central_region_width / 2)]

    if not gene_data.shape[0]:
        logging.warning(
            f"No data points available for slope fitting after positional fitering. Skipping dataset={dataset_id}")
        return [None, [np.nan, np.nan]]

    # slope_nc13 = np.nan
    # slope_nc14 = np.nan
    # max_polymerase_nc13, max_polymerase_nc14 = [np.nan] * 2
    # slope_start_pols = []

    # %% Calculate slopes and max values
    # plt.rc('text', usetex=False)
    set_figure_size(num=5, rows=1, page_width_frac=0.5, clear=True, height_factor=height_factor)
    fig, ax = plt.subplots(1, 1, num=5)
    # avg trace for a dataset

    avg_data, std_data = get_avg_on_regular_time_mesh(gene_data, dt=dt_new)
    x = avg_data.index
    y = avg_data.intensity
    ax.plot(x.values, y.values, '-o', markersize=markersize, lw=lw, c=colors['trace'])

    # Add error region
    ylow = y - std_data.intensity
    yup = y + std_data.intensity
    # print(ylow, yup)

    plt.fill_between(x.values, ylow.values, yup.values,
                     facecolor=colors['trace'], alpha=alpha, edgecolor=None)

    # x_text =
    # plt.text( f'nc {nc}')

    # ax.errorbar(x=avg_data.index.values, y=avg_data.intensity.values,
    #             yerr=std_data.values, fmt='-o', markersize=markersize, lw=lw)
    # print(std_data)
    # exit()

    # Detect ylims
    t_span = avg_data[~pd.isna(avg_data.intensity)].index.values
    xlims = np.array(t_span[[0, -1]])
    xlims = xlims + (xlims[1] - xlims[0]) * 0.1 / 2 * np.array([-1, 1])

    for nc in range(11, 15):
        nc_data = gene_data[gene_data.nc == nc]
        # print(nc_data)
        # return analyses
        if not nc_data.empty:
            # Add columns
            cols = ['slope', 'slopeV', 'max', 'maxV']
            for col in cols:
                if col not in analyses:
                    analyses[col] = np.nan

            start_slope = nc_data.time.min()
            end_slope = start_slope + slope_length_mins

            nc_data = gene_data[gene_data.nc == nc]

            # Calculate max and slope only if the nc includes at least 3 frames
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

            analyses.loc[(dataset_id, nc), cols] = (
                slope, slope_V, max_intensity, max_intensity_std**2)

            plot_slopes(gene_data, coefs, start_slope, end_slope, max_intensity, nc, ax)
    # print(analyses)
    plt.xlim(xlims)
    # save figure
    # print(dataset_name)
    figname = "slopes_{id:02d}_{dataset_name}".format(id=dataset_id, dataset_name=dataset_name)
    figpath = os.path.join(output_slopes_folder, figname)
    # print(figpath)

    if save_figures:
        # print('figpath', figpath)
        if pdf:
            fig.savefig(figpath + '.pdf', pad_inches=0, bbox_inches='tight')
        factor = 3
        figsize = fig.get_size_inches()
        fig.set_figwidth(figsize[0] * factor)
        fig.set_figheight(figsize[1] * factor)
        fig.savefig(figpath + '.png', pad_inches=0, bbox_inches='tight')

    # slope_series = pd.Series({'dataset_id': dataset_id, 'gene_id':
    #                           gene_id, 'construct_id': construct_id, 'slope_nc13': slope_nc13, 'slope_nc14': slope_nc14, 'max_nc13': max_polymerase_nc13, 'max_nc14': max_polymerase_nc14})  # , 'AP': AP})

    # print(f'Polymerase counts at the beginning of the slopes: {slope_start_pols.sort()}')
    # print(slope_start_pols)
    # print(slope_series)
    # print('Ola!')
    # %% Get a simple histogram of recorded AP positions
    fig = plt.figure(num=10, clear=True)
    # print('dat_name', dataset_name)
    str_title = '%s, id=%i' % (dataset_name.replace('_', '\_'), dataset_id)
    # print('str_tit', str_title)
    # if gene != "sn":

    # print(data.MeanAP.values)
    plt.hist(MeanAPs)
    plt.xlim([0, 1])

    plt.xlabel('AP')
    plt.title(str_title)

    figpath = f"AP_hist_{dataset_id}.png"
    figpath = os.path.join(AP_hist_folder, figpath)
    if save_figures:
        fig.savefig(figpath)
    # else:
    #     yPosMeans = data.groupby(data['trace_id'])[['yPos', 'yPosMean']].first()
    #     # print(yPosMeans.first())
    #
    #     plt.hist(yPosMeans.yPosMean, bins=range(0, 550, 50))
    #     # plt.xlim([0, 1])
    #
    #     plt.xlabel('yPosMean')
    #     plt.title(str_title)
    #
    #     figpath = f"AP_hist_{dataset_id}.png"
    #     figpath = os.path.join(AP_hist_folder, figpath)
    #     if save_figures:
    #         fig.savefig(figpath)

    return analyses


def get_slope(nc_data, start, end):
    cur_data = nc_data[(nc_data.time >= start)
                       & (nc_data.time <= end)]
    # print(nc, start, end)
    # print(cur_data)
    # print(len(cur_data))
    # print(cur_data.time.values, cur_data.intensity.values)
    if len(cur_data) <= 3:
        slope = np.nan
        coefs = [np.nan, np.nan]
        slope_V = np.nan
    else:
        coefs, V = np.polyfit(cur_data.time, cur_data.intensity, 1, cov=True)
        # print(slope)
        # print(V)
        slope = coefs[0]
        slope_V = V[0, 0]
        # print(slope_V)
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

    interval = max_time + np.array([-1, 1]) * dt_new
    int = nc_data[(nc_data.time >= max_time - dt_new) &
                  (nc_data.time <= max_time + dt_new)].intensity
    max_intensity = int.mean()
    max_intensity_std = np.std(int, ddof=1)
    return max_intensity, max_intensity_std


def plot_slopes(gene_data, coefs, start, end, max_intensity, nc, ax):
    dataset_name = gene_data.dataset.iloc[0]
    dataset_id = gene_data.dataset_id.iloc[0]

    # Plot slope fits
    # plot_data = gene_data.groupby(by="Frame").mean()
    # x_fit13 = np.asarray([start_nc13_frame, start_nc13_frame +
    #                       slope_length_frames - 1]) * mins_per_frame
    # y_fit13 = x_fit13 * coefs13[0] + coefs13[1]

    x_fit = np.asarray([start, end])
    y_fit = x_fit * coefs[0] + coefs[1]

    # # fit 13, 14
    # ax.plot(x_fit13, y_fit13, 'r-')
    ax.plot(x_fit, y_fit, 'r-', lw=lw, c=colors['slope'])

    # Expression threshold
    xlims = plt.xlim()
    intensity_threshold = intensity_thresholds.get(dataset_id, default_intensity_threshold)
    ax.plot(xlims, [intensity_threshold] * 2, 'g', lw=lw, c=colors['threshold'])

    # print nc number
    _, ymax = plt.ylim()
    plt.text(start, ymax * 0.9, f'nc{nc}', fontsize=7)

    #
    # # max values 13
    # end_nc13_frame = ncs_locations.iloc[dataset_id].nc13_end
    # x_max13 = np.asarray([start_nc13_frame, end_nc13_frame]) * mins_per_frame
    # y_max13 = np.asarray([1, 1]) * max_polymerase_nc13
    # plt.plot(x_max13, y_max13, 'g--')
    #
    # max values 14
    # start_nc14_frame = ncs_locations.iloc[dataset_id].nc14_start
    # end_nc14_frame = gene_data.Frame.max()
    x_max = gene_data[gene_data.nc == nc].time.min(), gene_data[gene_data.nc == nc].time.max()
    y_max = np.asarray([1, 1]) * max_intensity
    plt.plot(x_max, y_max, 'm--', lw=lw, c=colors['max'])

    # str_title = '%s, id = %i, dt = %.2f min' % (dataset_name, dataset_id, dt_new)
    # print('dat_name', dataset_name)
    str_title = dataset_name.replace('_', '\_')
    # print('str_title', str_title)
    plt.title(str_title)
    plt.xlabel('$t$, min')
    plt.ylabel('Fluo, a.u.')
    plt.ylim(ymin=-50)

    # plt.ylim([0, plt.ylim()[1]])


# def get_slope_old(gene_data, start_frame):
#     end_slope_frame = start_frame + slope_length_frames - 1
#     cur_data = gene_data[(gene_data.Frame >= start_frame) &
#                          (gene_data.Frame <= end_slope_frame)]
#     if cur_data.count().Frame <= 1 or np.isnan(start_frame):
#         slope = np.nan
#         coefs = [np.nan, np.nan]
#     else:
#         # print(cur_data)
#         coefs = np.polyfit(cur_data.time_min, cur_data.polymerases, 1)
#         slope = coefs[0]
#
#     # If the polymerase number at the start of the slope is too large, then we have missed the slope
#     slope_start_pol = cur_data[cur_data.Frame == start_frame].mean().polymerases
#
#     if slope_start_pol > 20:
#         slope_start_pol = np.nan
#         slope = np.nan
#         coefs = [np.nan, np.nan]
#
#     return [slope, coefs, slope_start_pol]


# def calculate_slopes_madhavs_data(data, ncs_locations):
#     """
#     Calculate nc13 and nc14 slope for the given dataset.
#     Calculate max polymerase numbers in these nc as well.
#     To minimize error, the slopes are calculated for the embryo as a whole from all points of the first frames.
#     This is more efficient than calculating for each trace independently (because it doesn't magnify the error), and more efficient than fitting just one straight line because it loses the frame wieghts.
#     """
#     dataset_name = data.dataset.iloc[0]
#     dataset_id = data.dataset_id.iloc[0]
#     gene_id = data.gene_id.iloc[0]
#     gene = data.gene.iloc[0]
#     construct_id = data.construct_id.iloc[0]
#     # AP = data.MeanAP.iloc[0]
#
#     # Filter out the highest gene-specific expression regions (one AP per trace)
#     MeanAPs = data.groupby(data['trace_id'])[['MeanAP']].first()
#     max_AP = MeanAPs.MeanAP.max()
#     median_AP = MeanAPs.MeanAP.median()
#
#     gene_data = data
#     if gene == 'hb':
#         # do not include the transition region of widht 0.05 at the right side
#         AP_trans_region_width = 0.05
#         gene_data = gene_data[(gene_data.MeanAP <= max_AP - AP_trans_region_width)
#                               & (gene_data.MeanAP >= 0.05)]
#     elif gene == 'kn':
#         AP_central_region_width = 0.04
#         gene_data = gene_data[(gene_data.MeanAP >= median_AP - AP_central_region_width / 2) &
#                               (gene_data.MeanAP <= median_AP + AP_central_region_width / 2)]
#
#     if not gene_data.shape[0]:
#         logging.warning(
#             f"No data points available for slope fitting after positional fitering. Skipping dataset={dataset_id}")
#         return [None, [np.nan, np.nan]]
#         # Mean AP
#     # AP = gene_data.MeanAP.mean()
#     # print(dataset_name)
#     # print([data.MeanAP.min(), data.MeanAP.max()])
#     # print(gene_data)
#     # print(AP)
#     # print(data)
#
#     # %% Get a simple histogram of recorded AP positions
#     fig = plt.figure(num=10, clear=True)
#     if gene != "sn":
#
#         # print(data.MeanAP.values)
#         plt.hist(MeanAPs.MeanAP)
#         plt.xlim([0, 1])
#
#         str_title = '%s, id=%i' % (dataset_name, dataset_id)
#
#         plt.xlabel('AP')
#         plt.title(str_title)
#
#         figpath = f"AP_hist_{dataset_id}.png"
#         figpath = os.path.join(AP_hist_folder, figpath)
#         fig.savefig(figpath)
#     else:
#         yPosMeans = data.groupby(data['trace_id'])[['yPos', 'yPosMean']].first()
#         # print(yPosMeans.first())
#
#         plt.hist(yPosMeans.yPosMean, bins=range(0, 550, 50))
#         # plt.xlim([0, 1])
#
#         str_title = '%s, id=%i' % (dataset_name, dataset_id)
#
#         plt.xlabel('yPosMean')
#         plt.title(str_title)
#
#         figpath = f"AP_hist_{dataset_id}.png"
#         figpath = os.path.join(AP_hist_folder, figpath)
#         fig.savefig(figpath)
#
#         # fig.savefig('AP_histogram.png')
#
#     # print(gene_data)
#
#     slope_start_pols = []
#     # nc 13 slope
#     start_nc13_frame = ncs_locations.iloc[dataset_id].nc13_start
#     # print(start_nc13_frame)
#     slope_nc13, coefs13, slope_start_pol = get_slope(gene_data, start_nc13_frame)
#     # print(slope_start_pol)
#
#     slope_start_pols.append(slope_start_pol)
#
#     # print(slope_nc13)
#
#     # nc 14 slope
#     start_nc14_frame = ncs_locations.iloc[dataset_id].nc14_start
#     slope_nc14, coefs14, slope_start_pol = get_slope(gene_data, start_nc14_frame)
#     slope_start_pols.append(slope_start_pol)
#     # print(slope_start_pols)
#
#     [max_polymerase_nc13, max_polymerase_nc14] = get_max_pol_number(
#         gene_data, dataset_id, ncs_locations)
#
#     # # Plot slope fits
#     # plot_data = gene_data.groupby(by="Frame").mean()
#     # x_fit13 = np.asarray([start_nc13_frame, start_nc13_frame +
#     #                       slope_length_frames - 1]) * mins_per_frame
#     # y_fit13 = x_fit13 * coefs13[0] + coefs13[1]
#     #
#     # x_fit14 = np.asarray([start_nc14_frame, start_nc14_frame +
#     #                       slope_length_frames - 1]) * mins_per_frame
#     # y_fit14 = x_fit14 * coefs14[0] + coefs14[1]
#     #
#     # plt.rc('text', usetex=False)
#     # fig = plt.figure(num=4)
#     # fig.clf()
#     # ax = fig.add_subplot(111)
#     #
#     # # data
#     # ax.plot(plot_data.time_min, plot_data.polymerases, '-o')
#     #
#     # # fit 13, 14
#     # ax.plot(x_fit13, y_fit13, 'r-')
#     # ax.plot(x_fit14, y_fit14, 'r-')
#     #
#     # # max values 13
#     # end_nc13_frame = ncs_locations.iloc[dataset_id].nc13_end
#     # x_max13 = np.asarray([start_nc13_frame, end_nc13_frame]) * mins_per_frame
#     # y_max13 = np.asarray([1, 1]) * max_polymerase_nc13
#     # plt.plot(x_max13, y_max13, 'g--')
#     #
#     # # max values 14
#     # start_nc14_frame = ncs_locations.iloc[dataset_id].nc14_start
#     # end_nc14_frame = gene_data.Frame.max()
#     # x_max = np.asarray([start_nc14_frame, end_nc14_frame]) * mins_per_frame
#     # y_max = np.asarray([1, 1]) * max_polymerase_nc14
#     # plt.plot(x_max, y_max, 'm--')
#     #
#     # str_title = '%s, id=%i' % (dataset_name, dataset_id)
#     # plt.title(str_title)
#     #
#     # plt.ylim([0, plt.ylim()[1]])
#     #
#     # # save figure
#     # figname = "slopes_%i.png" % dataset_id
#     # figpath = os.path.join(output_slopes_folder, figname)
#     # # print(figpath)
#     # fig.savefig(figpath)
#
#     slope_series = pd.Series({'dataset_id': dataset_id, 'gene_id':
#                               gene_id, 'construct_id': construct_id, 'slope_nc13': slope_nc13, 'slope_nc14': slope_nc14, 'max_nc13': max_polymerase_nc13, 'max_nc14': max_polymerase_nc14})  # , 'AP': AP})
#
#     # print(f'Polymerase counts at the beginning of the slopes: {slope_start_pols.sort()}')
#     # print(slope_start_pols)
#     # print(slope_series)
#     # print('Ola!')
#     return [slope_series, slope_start_pols]

# x = [-20.38971125, - 20.46165783, - 20.48963729, -
#      20.31975561]
# y = [387.51193237, 691.76159668, 376.16918945, 655.00183105]
# np.polyfit(x, y, 1, cov=True)

# print(np.__version__, sys.version)
