

import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# from constants import output_slopes_folder


def save_series_plot(x, y, a, b, fit_interval, vbars, filename, dataset_name, dataset_id):
    color_list = ['mediumorchid', 'yellowgreen', 'orange']

    fig = plt.figure(1)
    fig.clf()
    ax = plt.gca()

    # Plot trace
    ax.plot(x, y, '-o')

    # Plot fit
    # print([a, b])
    x_fit = np.asarray(fit_interval)
    y_fit = x_fit * a + b
    ax.plot(x_fit, y_fit, 'r')

    # Plot vertical bars
    y_bar = [np.min(y), np.max(y)]
    for i in range(len(vbars)):
        vbar = vbars[i]
        x_bar = [vbar, vbar]
        ax.plot(x_bar, y_bar, color=color_list[i])

    str_title = '%s, id=%i' % (dataset_name, dataset_id)
    # Calculate length
    nc13_length = vbars[1] - vbars[0]
    delta_nc = vbars[2] - vbars[0]
    if not np.isnan(nc13_length):
        str_title += ', nc13 = %.1f min' % nc13_length
    if not np.isnan(delta_nc):
        str_title += ', dnc = %.1f min' % delta_nc

    # Labels
    plt.xlabel("Time, min")
    plt.ylabel("Polymerase number")
    plt.title(str_title)

    # plt.show

    # Save to file
    # savepath = os.path.join(output_slopes_folder, filename)
    plt.savefig(filename)
