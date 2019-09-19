
import itertools
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constants import (LV, Darzacq2007, L, Tantale2016, colors_additivity,
                       figures_folder, k, kV, l, lV)
from least_squares_BCES import least_squares_BCES
from plot_theory_current_density import plot_theoretical_curve
from reinit_folder import reinit_folder
from set_figure_size import set_figure_size

height_factor = 0.65



def plot_normalized_current_density_diagram(analyses_in, I_est, num, I_estV=0, pdf=True):

    analyses = analyses_in.copy()
    long_labels = ['bac', 'no pr', 'no sh']
    gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}
    # colors = [[0, 0.251, 0], [0, 0, 1], [0.9412, 0.4706, 0], [0.502, 0.251, 0]]
    # ('b', 'brown', 'r', 'g'))

    genes = set(analyses.gene)
    constructs = set(analyses.construct)
    markers = ('x', 'o', '^', 'd', '*')
    # colors = ('#F16822', '#68972F', '#BE2026', 'mediumseagreen')
    colors = [colors_additivity[key] for key in colors_additivity]
    ms1 = 8
    ms2 = 3
    alpha = 0.25
    elinewidth = 1
    ylims = [0, 1.01]
    xlims = [0, 1.02]


    iterables = pd.MultiIndex.from_product([genes, constructs])
    Ts = pd.DataFrame(columns=['T', 'TV', 'Tstd', 'L', 'LV', 'Lstd'], index=iterables)

    # %% Plot a diagram with literature data
    fig = set_figure_size(num=num, rows=1, page_width_frac=0.5, height_factor=height_factor)
    ax = fig.subplots(1, 1)
    ax.clear()
    marker = itertools.cycle(markers)
    color = itertools.cycle(colors)

    # Plot
    for data in [Darzacq2007, Tantale2016]:
        c = 'k'  # next(color)
        m = next(marker)
        rho, rhoV, j, jV = [data[key] for key in ['rho', 'rhoV', 'j', 'jV']]
        ax.scatter(rho, j, s=ms1, c=c, marker=m,  zorder=4, label=data['label'])
        ax.errorbar(rho, j, xerr=rhoV**(1 / 2), yerr=jV**(1 / 2), fmt='.', markersize=ms2, c=c,
                    elinewidth=elinewidth, alpha=alpha, zorder=3)
    # plt.ylim(ymin=-1e-3)
    plot_theoretical_curve(ax, abort=False, norm=True)

    # Adjust & save
    def plot_adjust(analyses_filtered=None):
        plt.xlabel('Site occupation density $\\rho$')
        plt.ylabel('Normalized flux $j$')
        # plt
        # plt.xlim([0, 1.02])
        plt.xlim(xlims)
        plt.ylim(ylims)
        # plt.yticks(np.arange(0, 0.018, step=0.005))

        # ax.plot(plt.xlim(), np.array(plt.xlim()) / l, 'k', label=f'T = {T:.2f} min')
        if analyses_filtered is not None:
            T, TV, L, LV = analyses_filtered[['T', 'TV', 'L', 'LV']].head(1).values[0]
            L_rnd, Lstd_rnd = np.round(np.array([L, np.sqrt(LV)]) / 1e3, decimals=2)

            plt.title(
                f'{gene_long[gene]}')  # , $T = {T:.2f} \pm {np.sqrt(TV):.2f}$ min, $L = {L_rnd:.1f} \pm {Lstd_rnd:.1f}$ kb')

            # Add estimates
            ax.text(
                0.57, 0.05, f'$T = {T:.2f} \pm {np.sqrt(TV):.2f}$ min\n$L = {L_rnd:.1f} \pm {Lstd_rnd:.1f}$ kb', transform=ax.transAxes, va='bottom', ha='left', fontsize=8, weight='normal', style='normal', family='sans-serif')
        plt.legend(loc='upper left', labelspacing=0.2, frameon=False, borderaxespad=0)
        # else:
        # plt.legend()
        plt.tight_layout()
        fig.show()
    plot_adjust()
    plt.ylim(ymin=-0.02)

    figname = os.path.join(figures_folder, 'current-density_literature')
    fig.savefig(figname + '.png', pad_inches=0, bbox_inches='tight')
    if pdf:
        fig.savefig(figname + '.pdf', pad_inches=0, bbox_inches='tight')

    # %% Plot a diagram for our data
    num += 1
    for gene_id, gene in enumerate(genes):

        marker = itertools.cycle(markers)
        color = itertools.cycle(colors)

        set_figure_size(num=num, rows=1, page_width_frac=0.5, height_factor=height_factor)
        fig, ax = plt.subplots(1, 1, num=num, clear=True)
        markersize = 10

        for construct_id, construct in enumerate(constructs):

            # if gene == 'hb' and construct == 'bac':
            #     figlist = ['_without_data', '_with_data']
            # else:
            #     figlist = ['_with_data']
            figname = '_with_data'
            # ax.clear()

            analyses_filtered = analyses[(analyses.gene == gene)
                                         & (analyses.construct == construct)]

            # if figname == '_with_data':
            # for nc in range(11, 15):

            x = r = analyses_filtered.loc[:, 'rho'].values
            rV = analyses_filtered.loc[:, 'rhoV'].values
            y = j = analyses_filtered.loc[:, 'j'].values
            jV = analyses_filtered.loc[:, 'jV'].values

            xstd = np.sqrt(rV)
            ystd = np.sqrt(jV)

            c = next(color)
            m = next(marker)
            ax.scatter(x, y,
                       label=f'{long_labels[construct_id]}', s=ms1, c=c, marker=m,  zorder=4)
            ax.errorbar(x, y, xerr=xstd, yerr=ystd, fmt='.', markersize=ms2, c=c,
                        elinewidth=elinewidth, alpha=alpha, zorder=3)

            plot_theoretical_curve(ax, abort=False, norm=True)
            plot_adjust(analyses_filtered)

            # ax.plot([gamma_max] * 2, ylims, '-r')
            # ax.plot([1] * 2, ylims, '-k')
            plt.ylim(ylims)
            plt.xlim(xlims)

            figname = os.path.join(figures_folder, 'current-density_' +
                                   gene_long[gene] + figname)
            fig.savefig(figname + '.png', pad_inches=0, bbox_inches='tight')
            if pdf:
                fig.savefig(figname + '.pdf', pad_inches=0, bbox_inches='tight')

    return analyses
