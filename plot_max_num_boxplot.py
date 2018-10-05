

import os

import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from constants import (construct_labels, gene_labels, hb_AP_interval,
                       kn_AP_interval, max_pol_abortive_theory,
                       max_pol_simple_theory, slope_abortive_theory,
                       slope_simple_theory)
from set_figure_size import set_figure_size


def plot_max_num_boxplot(slopes):

    height_factor = 3
    markersize = 2
    ylims = [0, 120]
    hist_color = np.asarray([117, 145, 41]) / 255
    # hist_color = "#89bc00"
    # color1 = "#9e00f8"
    # color2 = "#f80000"
    color1 = 'k'
    color2 = 'k'
    alpha = 0.1
    label_y = 0.925
    font_size = 8
    dashes = (6, 6)
    linewidth = 0.5
    box_width = 0.1
    xlims = [0.5, 18.5]
    counts_y = 25

    fig = set_figure_size(num=3, rows=1, page_width_frac=0.5, height_factor=height_factor)

    nc_labels = ['max_nc13', 'max_nc14']

    gene_names_full = ['Hunchback', 'Knirps', 'Snail']

    datasets = []
    datasets_labels = []
    bp = None
    counts = []

    for gene_id in range(3):
        slopes_gene = slopes[(slopes.gene_id == gene_id)]
        for nc in [13, 14]:
            for construct_id in [0, 2, 1]:
                nc_label = nc_labels[nc - 13]
                dataset = slopes_gene[(
                    slopes_gene.construct_id == construct_id)][nc_label].dropna()
                count = dataset.count()

                # if gene_id == 1 and construct_id in [2, 1]:
                #     print(nc, construct_id)
                #     print(dataset)

                datasets.append(dataset.values)
                datasets_labels.append(f"{construct_labels[construct_id]}{nc}")
                counts.append(count)

    # Plot
    ax = fig.add_subplot(111)

    # Add gene rectangles and labels
    # Hb
    rect = patches.Rectangle((0, ylims[0]), 18 / 3 + 0.5, ylims[1] - ylims[0],
                             edgecolor='none', facecolor=color1, alpha=alpha)
    ax.add_patch(rect)
    ax.text(1 / 6, label_y, gene_names_full[0],
            transform=ax.transAxes, fontsize=font_size, ha='center')

    # Kn
    ax.text(1 / 2, label_y, gene_names_full[1],
            transform=ax.transAxes, fontsize=font_size, ha='center')

    # Sn
    rect = patches.Rectangle((0.5 + 18 * 2 / 3, ylims[0]), 18, ylims[1] - ylims[0],
                             edgecolor='none', facecolor=color2, alpha=alpha)
    ax.add_patch(rect)
    ax.text(1 - 1 / 6, label_y, gene_names_full[2],
            transform=ax.transAxes, fontsize=font_size, ha='center')

    # # Add number of datasets
    # for i in range(18):
    #     ax.text(i + 1, counts_y, counts[i], fontsize=font_size, ha='center')

    # # plot the scatter plot of the data itself
    # for i in range(len(datasets)):
    #     y = datasets[i]
    #     x = np.random.normal(i + 1, box_width, size=len(y))
    #     ax.plot(x, y, 'ok', markersize=markersize, alpha=0.3)

    # Plot the boxplot
    bp = ax.boxplot(datasets, showfliers=True)

    # Adjust
    plt.ylim(ylims)
    plt.xlim(xlims)

    # Labels
    ax.set_xticks(range(1, 19))
    ax.set_xticklabels(datasets_labels)
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right",
             rotation_mode="anchor")

    # change the style of fliers
    if bp != None:
        for flier in bp['fliers']:
            flier.set(marker='o', color='#e7298a', alpha=0.5, markersize=markersize)
            # print(flier)

    # Add theory
    x_theor = plt.xlim()
    y_theor_simple = np.asarray([1, 1]) * max_pol_simple_theory
    y_theor_abortive = np.asarray([1, 1]) * max_pol_abortive_theory
    print(max_pol_abortive_theory)

    ax.plot(x_theor, y_theor_simple, 'r-', linewidth=linewidth)  # , dashes=dashes)
    ax.plot(x_theor, y_theor_abortive, 'b--', linewidth=linewidth, dashes=dashes)

    # Labels
    plt.ylabel('Max. polymerase number')

    plt.show()
