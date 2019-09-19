
import itertools
import os

import matplotlib.pyplot as plt
import numpy as np

from constants import (alpha_additivity, colors_additivity, figures_folder,
                       markers_additivity)
from set_figure_size import set_figure_size
from theoretical_phase_parameters import (J_over_k_LD, alpha_over_k_abortive,
                                          rho_LD)


def plot_parameter_evolution(analyses, gene='hb', pdf=False):

    # all_constructs = ['bac', 'no_pr', 'no_sh']
    constructs = set(analyses.construct)

    # long_labels = ['no shadow', 'no primary', 'bac'] A2142F
    long_labels = {'bac': 'bac', 'no_pr': 'no pr', 'no_sh': 'no sh'}
    gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}
    y_label = {'j': 'Normalized flux $j$',
               'rho': 'Site occupation density $\\rho$', 'tau': 'Residence time $\\tau$ (s)', 'alpha_comb': 'Initiation rate $\\alpha$ (pol/min)'}
    x_jiggle = 0.04
    x_shifts = np.array([-1, 0, 1]) * x_jiggle
    capsize = 0
    markersize = 4
    lw = 1
    ncs = np.arange(11, 15)
    p_shift = 0.12

    grouped_data = analyses.groupby(by=['gene', 'construct', 'nc'])
    all_means = grouped_data.mean()
    all_stds = grouped_data.std(ddof=1)
    all_ps = analyses.groupby(by=['gene', 'nc']).first()
    # print(all_ps)

    for j, quantity in enumerate(['j', 'rho', 'tau', 'alpha_comb']):
        # print(quantity)
        ymaxs = {'j': 0.36, 'rho': 0.27, 'tau': 103, 'alpha_comb': 12}
        num = 12 + j
        set_figure_size(num=num, rows=1, page_width_frac=0.5, clear=True, height_factor=0.7)
        fig, ax = plt.subplots(1, 1, num=num, clear=True)
        # label = f'{quantity}_comb'
        avg_data, std_data = {}, {}
        for construct in constructs:
            if quantity in ['rho', 'j']:
                avg_data[construct] = all_means.loc[(
                    gene, construct, slice(None)), quantity].values

                std_data[construct] = all_stds.loc[(gene, construct, slice(None)), quantity].values
            elif quantity in ['tau', 'alpha_comb']:
                avg_data[construct] = all_means.loc[(
                    gene, construct, slice(None)), quantity].values
                std_data[construct] = np.sqrt(
                    all_means.loc[(gene, construct, slice(None)), quantity + 'V'].values)

        marker_gen = itertools.cycle(markers_additivity)

        for i, construct in enumerate(constructs):
            # print(colors_additivity[construct])
            m = next(marker_gen)
            plt.errorbar(ncs + x_shifts[i], avg_data[construct], yerr=std_data[construct],
                         fmt='-' + m,  color=colors_additivity[construct], capsize=capsize, label=long_labels[construct], markersize=markersize, lw=lw)

        # if gene == 'kn':
        #     print(quantity, avg_data, std_data)

        # Make a filled region for the sum
        x_fill = ncs.copy().astype(np.float)
        # print(x_fill)
        x_fill[[0, 1]] -= x_jiggle
        x_fill[[-2, -1]] += x_jiggle

        plt.xlabel('Nuclear cycle')
        plt.ylabel(y_label[quantity])
        plt.ylim(ymin=0)

        plt.ylim(ymax=ymaxs[quantity])

        # Add p-values
        if quantity in []:  # ['alpha_comb']:
            p = all_ps.loc[(gene, slice(None)), quantity + '_p_value'].values
            # y = np.nanmax(avg_data[construct] + std_data[construct]) * 1.1
            # ymax = y * 1.4
            ymax = ymaxs[quantity]
            # if quantity is 'tau':
            # ymax = y / 1.05
            y = ymax / 1.4
            # if quantity is not 'tau':
            for i, nc in enumerate(ncs):
                x = nc - p_shift
                # if quantity == 'JoK':
                #     y = 6e-3
                # elif quantity == 'rho':
                #     y = 0.27

                # y = abrt_limit * 1.2

                # p_rnd =
                if not np.isnan(p[i]):
                    digits = -int(np.floor(np.log10(p[i])))
                    p_rnd = np.round(p[i], decimals=np.max([2, digits]))
                    ax.text(x, y, f'p={p_rnd}', rotation=90, va='bottom', ha='center')
                # plt.ylim(ymax=ymax)
        # print(p)

        # _, ymax = plt.ylim()
        # if quantity not in ['tau']:

        # print(ymax)

        # else:
        # plt.ylim(ymax=60)
        # else:

        plt.xticks(ncs)
        plt.title(gene_long[gene])
        # if quantity == 'tau':
        #     plt.legend(loc='upper right', labelspacing=0.2, frameon=False, borderaxespad=0)
        plt.tight_layout()
        plt.show()

        figname = 'additivity_' + quantity + '_' + gene
        figpath = os.path.join(figures_folder, figname)
        fig.savefig(figpath + '.png', pad_inches=0, bbox_inches='tight')
        if pdf:
            fig.savefig(figpath + '.pdf', pad_inches=0, bbox_inches='tight')
