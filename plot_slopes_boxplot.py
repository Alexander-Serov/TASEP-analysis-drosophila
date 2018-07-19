

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
from set_figure_size import set_figure_size

from constants import slope_simple_theory, hb_AP_interval, kn_AP_interval


def plot_slopes_boxplot(slopes):

    height_factor = 3
    markersize = 2
    ylims = [-5, 30]
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

    fig = set_figure_size(num=2, rows=1, page_width_frac=0.5, height_factor=height_factor)

    nc_labels = ['slope_nc13', 'slope_nc14']
    construct_labels = ['b', 's', 'p']
    genes_labels = ['Hunchback', 'Knirps', 'Snail']

    datasets = []
    datasets_labels = []
    bp = None

    for gene_id in range(3):
        # Keep only maximal expression regions for each gene
        slopes_gene = slopes[(slopes.gene_id == gene_id)]
        if gene_id == 0:  # hb
            slopes_gene = slopes[(slopes.gene_id == gene_id) & (
                slopes.AP >= hb_AP_interval[0]) & (slopes.AP <= hb_AP_interval[1])]
        elif gene_id == 1:  # kn
            slopes_gene = slopes[(slopes.gene_id == gene_id) & (
                slopes.AP >= kn_AP_interval[0]) & (slopes.AP <= kn_AP_interval[1])]

        for nc in [13, 14]:
            for construct_id in [0, 2, 1]:

                nc_label = nc_labels[nc - 13]

                # Average the data by the embryo
                group_data = slopes_gene[(
                    slopes_gene.construct_id == construct_id)].groupby(by="dataset_id")
                # print(group_data
                group_data = group_data[nc_label].mean()
                group_data = group_data.dropna().values

                print(len(group_data))
                datasets.append(group_data)
                datasets_labels.append("%s%i" % (construct_labels[construct_id], nc))

    # print(datasets)

    # Plot
    ax = fig.add_subplot(111)

    # Add gene rectangles and labels
    # Hb
    rect = patches.Rectangle((0, ylims[0]), 18 / 3 + 0.5, ylims[1] - ylims[0],
                             edgecolor='none', facecolor=color1, alpha=alpha)
    ax.add_patch(rect)
    ax.text(1 / 6, label_y, genes_labels[0],
            transform=ax.transAxes, fontsize=font_size, ha='center')

    # Kn
    ax.text(1 / 2, label_y, genes_labels[1],
            transform=ax.transAxes, fontsize=font_size, ha='center')

    # Sn
    rect = patches.Rectangle((0.5 + 18 * 2 / 3, ylims[0]), 18, ylims[1] - ylims[0],
                             edgecolor='none', facecolor=color2, alpha=alpha)
    ax.add_patch(rect)
    ax.text(1 - 1 / 6, label_y, genes_labels[2],
            transform=ax.transAxes, fontsize=font_size, ha='center')

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
    y_theor = np.asarray([1, 1]) * slope_simple_theory
    ax.plot(x_theor, y_theor, 'k--', linewidth=linewidth, dashes=dashes)

    plt.show()
