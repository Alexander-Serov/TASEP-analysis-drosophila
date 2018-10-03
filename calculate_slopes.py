

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from constants import slope_length_frames, hb_AP_interval, kn_AP_interval, mins_per_frame, output_slopes_folder


def calculate_slopes(data, ncs_locations):
    """
    Calculate nc13 and nc14 slope for the given dataset
    """
    dataset_name = data.dataset.iloc[0]
    dataset_id = data.dataset_id.iloc[0]
    gene_id = data.gene_id.iloc[0]
    construct_id = data.construct_id.iloc[0]
    # AP = data.MeanAP.iloc[0]

    # Filter out the highest gene-specific expression regions
    gene_data = data
    if gene_id == 0:
        gene_data = gene_data[(gene_data.MeanAP >= hb_AP_interval[0]) &
                              (gene_data.MeanAP <= hb_AP_interval[1])]
    elif gene_id == 1:
        gene_data = gene_data[(gene_data.MeanAP >= kn_AP_interval[0]) &
                              (gene_data.MeanAP <= kn_AP_interval[1])]

    # Mean AP
    AP = gene_data.MeanAP.mean()

    def get_slope(start_frame):
        end_slope_frame = start_frame + slope_length_frames - 1
        cur_data = gene_data[(gene_data.Frame >= start_frame) &
                             (gene_data.Frame <= end_slope_frame)]
        if cur_data.count().Frame <= 1 or np.isnan(start_frame):
            slope = np.nan
            coefs = [np.nan, np.nan]
        else:
            coefs = np.polyfit(cur_data.time_min, cur_data.polymerases, 1)
            slope = coefs[0]

        # Discard slope if negative
        if slope < 0:
            slope = np.nan
        return [slope, coefs]

    # print(gene_data)

    # nc 13 slope
    start_nc13_frame = ncs_locations.iloc[dataset_id].nc13_start
    # print(start_nc13_frame)
    slope_nc13, coefs13 = get_slope(start_nc13_frame)
    # print(slope_nc13)

    # nc 14 slope
    start_nc14_frame = ncs_locations.iloc[dataset_id].nc14_start
    slope_nc14, coefs14 = get_slope(start_nc14_frame)

    max_polymerase_nc13 = np.nan
    max_polymerase_nc14 = np.nan

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

    str_title = '%s, id=%i' % (dataset_name, dataset_id)
    plt.title(str_title)

    # save figure
    figname = "slopes_dataset_%i.png" % dataset_id
    figpath = os.path.join(output_slopes_folder, figname)
    fig.savefig(figpath)

    # plt.show()

    # dataset_id

    #
    # def get_mean_max(start_frame, end_frame):
    #     cur_data = data[(data.Frame >= start_frame) & (data.Frame <= end_frame)]
    #     if cur_data.count().Frame < 3 or np.isnan(start_frame) or np.isnan(end_frame):
    #         max_polymerase = np.nan
    #     else:
    #         max_frame = cur_data[cur_data.polymerases ==
    #                              cur_data.polymerases.max()].Frame.iloc[0]
    #         # Get the average in 3 frames around the maximum. If there are no frames around, skip
    #         data_around_max = cur_data[(cur_data.Frame >= max_frame - 1) &
    #                                    (cur_data.Frame <= max_frame + 1)]
    #         if data_around_max.count().Frame < 3:
    #             max_polymerase = np.nan
    #         else:
    #             max_polymerase = data_around_max.mean().polymerases
    #     return max_polymerase

    # nc 13 slope

    #
    # # nc 13 max value
    # end_nc13_frame = ncs_locations.iloc[dataset_id].nc13_end
    # max_polymerase_nc13 = get_mean_max(start_nc13_frame, end_nc13_frame)
    #

    #
    # # nc 14 max value
    # end_nc14_frame = data[(data.trace_id == trace_id)].max().Frame
    # max_polymerase_nc14 = get_mean_max(start_nc14_frame, end_nc14_frame)
    #
    slope_series = pd.Series({'dataset_id': dataset_id, 'gene_id':
                              gene_id, 'construct_id': construct_id, 'slope_nc13': slope_nc13, 'slope_nc14': slope_nc14, 'max_nc13': max_polymerase_nc13, 'max_nc14': max_polymerase_nc14, 'AP': AP})
    return slope_series
