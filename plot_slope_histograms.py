

import matplotlib.pyplot as plt
import numpy as np
from set_figure_size import set_figure_size

from constants import slope_theory


def plot_slope_histograms(slopes):
    bin_width = 1.75
    rows = 1
    cols = 3
    alpha = 0.8
    xlims = [-15, 30]
    genes_names = ['Hunchback', 'Knirps', 'Snail']
    hist_color = np.asarray([117, 145, 41]) / 255
    linewidth = 0.5
    dashes = (6, 6)

    # fig = plt.figure(num=2)
    # fig.clf()
    set_figure_size(num=2, rows=1, page_width_frac=1, height_factor=0.8)
    fig, axarr = plt.subplots(num=2, nrows=rows, ncols=cols, sharex=True, sharey=True)

    def plot_one_hist(slopes_no_nan, **kwargs):
        # slopes_no_nan = cur_slopes.dropna()
        slopes_no_nan.min()
        bins = np.arange(slopes_no_nan.min(), slopes_no_nan.max() + bin_width, bin_width)
        plt.hist(slopes_no_nan, bins=bins, density=True, histtype='bar', ec='black', linewidth=linewidth,
                 alpha=alpha, color=hist_color, ** kwargs)

    for gene_id in range(3):
        # for construct_id in range(3):

        ax = axarr[gene_id]
        plt.sca(ax)

        gene_slopes = slopes[(slopes.gene_id == gene_id)]
        gene_slopes_all_nc = np.concatenate(
            [gene_slopes.slope_nc13.dropna().values, gene_slopes.slope_nc14.dropna().values])
        traces_len = len(gene_slopes_all_nc)
        plot_one_hist(gene_slopes_all_nc)

        # Adjust
        plt.xlabel('Initial injection rate (pol/min)')
        if gene_id == 0:
            plt.ylabel('PDF')
        plt.xlim(xlims)
        str_title = '\emph{%s} ($n=%i$)' % (genes_names[gene_id], traces_len)
        plt.title(str_title)

    # Add theoretical limit
    x_theor = np.asarray([1, 1]) * slope_theory
    y_theor = np.asarray(plt.ylim())
    # y_theor[1] *= 0.9
    for ax in axarr:
        ax.plot(x_theor, y_theor, 'k--', linewidth=linewidth, dashes=dashes)

    fig.tight_layout()
    plt.show()
