

import matplotlib.pyplot as plt
import numpy as np

from constants import (gene_labels, max_pol_abortive_theory,
                       max_pol_simple_theory, slope_abortive_theory,
                       slope_simple_theory)
from set_figure_size import set_figure_size


def plot_boxplot_hb(slopes):

    height_factor = 1
    markersize = 2.5
    # ylims = [0, 120]
    # hist_color = np.asarray([117, 145, 41]) / 255
    # color1 = 'k'
    # color2 = 'k'
    # alpha = 0.1
    # label_y = 0.925
    # font_size = 8
    dashes = (6, 6)
    linewidth = 0.5
    # box_width = 0.1
    # xlims = [0.5, 18.5]
    # counts_y = 25
    gene_id = 0
    # ylims_slope = [0, 30]
    # ylims_number = [0, 120]

    fig = set_figure_size(num=3, rows=1, page_width_frac=1, height_factor=height_factor)

    nc_labels = ['max_nc13', 'max_nc14']
    y_labels = {'slope': 'Initial injection rate (pol / min)',
                'max': 'Max. polymerase number'}
    ylims = {'slope': [0, 30], 'max': [0, 120]}
    theor_simple = {'slope': slope_simple_theory, 'max': max_pol_simple_theory}
    theor_abortive = {'slope': slope_abortive_theory, 'max': max_pol_abortive_theory}

    datasets = []
    datasets_labels = ['bac', 'no sh', 'no pr']
    # bp = None
    # counts = []

    gene_data = slopes[(slopes.gene_id == gene_id)]
    for nc in [13, 14]:
        for construct_id in [0, 2, 1]:
            nc_label = nc_labels[nc - 13]
            dataset = gene_data[(
                gene_data.construct_id == construct_id)][nc_label].dropna()

            datasets.append(dataset.values)
            # datasets_labels.append(f"{construct_labels[construct_id]}{nc}")

    # %% Plot
    fig.clf()
    axes = fig.subplots(nrows=1, ncols=4)
    axes_gen = (axis for axis in axes)

    def make_one_plot(type_, nc):
        label = type_ + '_nc' + str(nc)
        ax = next(axes_gen)
        datasets = []
        for construct_id in [0, 2, 1]:
            cur_data = gene_data[(gene_data.construct_id == construct_id)][label].dropna()
            datasets.append(cur_data)
            # datasets_labels.append(f"{construct_labels[construct_id]}{nc}")
        bp = ax.boxplot(datasets, showfliers=True)

        # Adjust
        # ax.set_ylim([0, ax.get_ylim()[1]])
        ax.set_ylim(ylims[type_])
        ax.set_xticklabels(datasets_labels)
        plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
                 rotation_mode="anchor")
        ax.set_title(f"nc{nc}")
        ax.set_ylabel(y_labels[type_])

        # change the style of fliers
        if bp is not None:
            for flier in bp['fliers']:
                flier.set(marker='o', color='#e7298a', alpha=0.5, markersize=markersize)

        # Add theory
        x_theor = ax.get_xlim()
        y_theor_simple = np.asarray([1, 1]) * theor_simple[type_]
        y_theor_abortive = np.asarray([1, 1]) * theor_abortive[type_]
        ax.plot(x_theor, y_theor_simple, 'r-', linewidth=linewidth)  # , dashes=dashes)
        ax.plot(x_theor, y_theor_abortive, 'b--', linewidth=linewidth, dashes=dashes)

        return ax

    # slopes
    make_one_plot('slope', 13)
    # ax.set_ylabel('Initial injection slope, pol / min')
    make_one_plot('slope', 14)
    # max. pol. numbers
    make_one_plot('max', 13)
    # ax.set_ylabel('Max. polymerase number')
    make_one_plot('max', 14)

    # plt.xlim(xlims)

    # Labels
    # ax.set_xticks(range(1, 19))

    # print(flier)

    # Add theory

    #
    # print(max_pol_abortive_theory)

    #

    # Labels
    # plt.ylabel('Max. polymerase number')
    # plt.title(f"nc{nc}")

    fig.tight_layout()

    plt.show()
