

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os


from constants import output_slopes_folder


def save_series_plot(x, y, a, b, fit_interval, filename, dataset_name):

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

    # Labels
    plt.xlabel("Time, min")
    plt.xlabel("Polymerase number")
    plt.title(dataset_name)

    # plt.show

    # Save to file
    savepath = os.path.join(output_slopes_folder, filename)
    plt.savefig(savepath)
