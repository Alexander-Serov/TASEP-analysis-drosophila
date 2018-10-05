

import logging
import os
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constants import (AP_hist_folder, hb_AP_interval, kn_AP_interval,
                       mins_per_frame, output_slopes_folder,
                       slope_length_frames)


def calculate_slopes(data, ncs_locations):
    """
    Calculate nc13 and nc14 slope for the given dataset.
    Calculate max polymerase numbers in these nc as well.
    To minimize error, the slopes are calculated for the embryo as a whole from all points of the first frames.
    This is more efficient than calculating for each trace independently (because it doesn't magnify the error), and more efficient than fitting just one straight line because it loses the frame wieghts.
    """
    dataset_name = data.dataset.iloc[0]
    dataset_id = data.dataset_id.iloc[0]
    gene_id = data.gene_id.iloc[0]
    gene = data.gene.iloc[0]
    construct_id = data.construct_id.iloc[0]
    # AP = data.MeanAP.iloc[0]

    # Filter out the highest gene-specific expression regions (one AP per trace)
    MeanAPs = data.groupby(data['trace_id'])[['MeanAP']].first()
    max_AP = MeanAPs.MeanAP.max()
    median_AP = MeanAPs.MeanAP.median()

    gene_data = data
    if gene == 'hb':
        # do not include the transition region of widht 0.05 at the right side
        AP_trans_region_width = 0.05
        gene_data = gene_data[(gene_data.MeanAP <= max_AP - AP_trans_region_width)
                              & (gene_data.MeanAP >= 0.05)]
    elif gene == 'kn':
        AP_central_region_width = 0.04
        gene_data = gene_data[(gene_data.MeanAP >= median_AP - AP_central_region_width / 2) &
                              (gene_data.MeanAP <= median_AP + AP_central_region_width / 2)]

    if not gene_data.shape[0]:
        logging.warning(
            f"No data points available for slope fitting after positional fitering. Skipping dataset={dataset_id}")
        return [None, [np.nan, np.nan]]
        # Mean AP
    # AP = gene_data.MeanAP.mean()
    # print(dataset_name)
    # print([data.MeanAP.min(), data.MeanAP.max()])
    # print(gene_data)
    # print(AP)
    # print(data)

    # %% Get a simple histogram of recorded AP positions
    fig = plt.figure(num=10, clear=True)
    if gene != "sn":

        # print(data.MeanAP.values)
        plt.hist(MeanAPs.MeanAP)
        plt.xlim([0, 1])

        str_title = '%s, id=%i' % (dataset_name, dataset_id)

        plt.xlabel('AP')
        plt.title(str_title)

        figpath = f"AP_hist_{dataset_id}.png"
        figpath = os.path.join(AP_hist_folder, figpath)
        fig.savefig(figpath)
    else:
        yPosMeans = data.groupby(data['trace_id'])[['yPos', 'yPosMean']].first()
        # print(yPosMeans.first())

        plt.hist(yPosMeans.yPosMean, bins=range(0, 550, 50))
        # plt.xlim([0, 1])

        str_title = '%s, id=%i' % (dataset_name, dataset_id)

        plt.xlabel('yPosMean')
        plt.title(str_title)

        figpath = f"AP_hist_{dataset_id}.png"
        figpath = os.path.join(AP_hist_folder, figpath)
        fig.savefig(figpath)

        # fig.savefig('AP_histogram.png')

    # print(gene_data)

    slope_start_pols = []
    # nc 13 slope
    start_nc13_frame = ncs_locations.iloc[dataset_id].nc13_start
    # print(start_nc13_frame)
    slope_nc13, coefs13, slope_start_pol = get_slope(gene_data, start_nc13_frame)
    # print(slope_start_pol)

    slope_start_pols.append(slope_start_pol)

    # print(slope_nc13)

    # nc 14 slope
    start_nc14_frame = ncs_locations.iloc[dataset_id].nc14_start
    slope_nc14, coefs14, slope_start_pol = get_slope(gene_data, start_nc14_frame)
    slope_start_pols.append(slope_start_pol)
    # print(slope_start_pols)

    [max_polymerase_nc13, max_polymerase_nc14] = get_max_pol_number(
        gene_data, dataset_id, ncs_locations)

    # Plot slope fits
    plot_data = gene_data.groupby(by="Frame").mean()
    x_fit13 = np.asarray([start_nc13_frame, start_nc13_frame +
                          slope_length_frames - 1]) * mins_per_frame
    y_fit13 = x_fit13 * coefs13[0] + coefs13[1]

    x_fit14 = np.asarray([start_nc14_frame, start_nc14_frame +
                          slope_length_frames - 1]) * mins_per_frame
    y_fit14 = x_fit14 * coefs14[0] + coefs14[1]

    plt.rc('text', usetex=False)
    fig = plt.figure(num=4)
    fig.clf()
    ax = fig.add_subplot(111)

    # data
    ax.plot(plot_data.time_min, plot_data.polymerases, '-o')

    # fit 13, 14
    ax.plot(x_fit13, y_fit13, 'r-')
    ax.plot(x_fit14, y_fit14, 'r-')

    # max values 13
    end_nc13_frame = ncs_locations.iloc[dataset_id].nc13_end
    x_max13 = np.asarray([start_nc13_frame, end_nc13_frame]) * mins_per_frame
    y_max13 = np.asarray([1, 1]) * max_polymerase_nc13
    plt.plot(x_max13, y_max13, 'g--')

    # max values 14
    start_nc14_frame = ncs_locations.iloc[dataset_id].nc14_start
    end_nc14_frame = gene_data.Frame.max()
    x_max = np.asarray([start_nc14_frame, end_nc14_frame]) * mins_per_frame
    y_max = np.asarray([1, 1]) * max_polymerase_nc14
    plt.plot(x_max, y_max, 'm--')

    str_title = '%s, id=%i' % (dataset_name, dataset_id)
    plt.title(str_title)

    plt.ylim([0, plt.ylim()[1]])

    # save figure
    figname = "slopes_%i.png" % dataset_id
    figpath = os.path.join(output_slopes_folder, figname)
    # print(figpath)
    fig.savefig(figpath)

    slope_series = pd.Series({'dataset_id': dataset_id, 'gene_id':
                              gene_id, 'construct_id': construct_id, 'slope_nc13': slope_nc13, 'slope_nc14': slope_nc14, 'max_nc13': max_polymerase_nc13, 'max_nc14': max_polymerase_nc14})  # , 'AP': AP})

    # print(f'Polymerase counts at the beginning of the slopes: {slope_start_pols.sort()}')
    # print(slope_start_pols)
    # print(slope_series)
    # print('Ola!')
    return [slope_series, slope_start_pols]


def get_slope(gene_data, start_frame):
    end_slope_frame = start_frame + slope_length_frames - 1
    cur_data = gene_data[(gene_data.Frame >= start_frame) &
                         (gene_data.Frame <= end_slope_frame)]
    if cur_data.count().Frame <= 1 or np.isnan(start_frame):
        slope = np.nan
        coefs = [np.nan, np.nan]
    else:
        # print(cur_data)
        coefs = np.polyfit(cur_data.time_min, cur_data.polymerases, 1)
        slope = coefs[0]

    # If the polymerase number at the start of the slope is too large, then we have missed the slope
    slope_start_pol = cur_data[cur_data.Frame == start_frame].mean().polymerases

    if slope_start_pol > 20:
        slope_start_pol = np.nan
        slope = np.nan
        coefs = [np.nan, np.nan]

    return [slope, coefs, slope_start_pol]


def get_max_pol_number(gene_data, dataset_id, ncs_locations):
    """
    Get the maximum polymerase number from the cycle. Average over +-1 frame for less noise.
    """
    dataset = gene_data[gene_data.dataset_id == dataset_id].dataset.iloc[0]
    max_polymerase_nc13, max_polymerase_nc14 = np.nan, np.nan
    # nc 13
    start_frame = ncs_locations.iloc[dataset_id].nc13_start
    end_frame = ncs_locations.iloc[dataset_id].nc13_end

    # manual corrections
    if dataset == '2014-06-08-SnaBACA':
        end_frame -= 20

    if np.all(np.isfinite([start_frame, end_frame])):
        nc_data = gene_data[(gene_data.Frame >= start_frame) &
                            (gene_data.Frame <= end_frame)]
        nc_data = nc_data.groupby(by='Frame').mean()
        # print(nc_data.polymerases == nc_data.polymerases.max())
        max_frame = nc_data.polymerases.idxmax()
        max_polymerase_nc13 = nc_data.loc[max_frame - 1:max_frame + 1].polymerases.mean()

    # nc 14
    start_frame = ncs_locations.iloc[dataset_id].nc14_start
    if np.isfinite(start_frame):
        nc_data = gene_data[(gene_data.Frame >= start_frame)]
        nc_data = nc_data.groupby(by='Frame').mean()
        max_frame = nc_data.polymerases.idxmax()
        max_polymerase_nc14 = nc_data.loc[max_frame - 1:max_frame + 1].polymerases.mean()

    # print([max_polymerase_nc13, max_polymerase_nc14])
    return [max_polymerase_nc13, max_polymerase_nc14]
