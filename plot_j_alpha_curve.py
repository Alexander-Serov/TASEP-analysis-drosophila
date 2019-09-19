import itertools
import os

import numpy as np
from matplotlib import pyplot as plt

from constants import (Darzacq2007, Tantale2016, colors_additivity,
                       figures_folder, gene_long, l)
from set_figure_size import set_figure_size

# Constants
alpha_over_k_min = 0
alpha_over_k_max = 1 / (1 + np.sqrt(l))  # theoretical maximum
steps = 100
ms1 = 8
ms2 = 4  # size of the marker for errorbar plot to superimpose on the the scatter data plot
alpha = 0.25
height_factor = 0.65
elinewidth = 1  # errorbar line width
markers = ('x', 'o', '^', 'd', '*')
colors_theor = [colors_additivity['sum'],  colors_additivity['brown']]
colors = [colors_additivity[key] for key in colors_additivity]
long_labels = ['bac', 'no pr', 'no sh']
gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}
genes = ['hb', 'sn', 'kn']
constructs = ['bac', 'no_pr', 'no_sh']
figname = 'j_alpha'


def plot_theory(ax):
    """Plot theoretical curves for rho and j"""
    alpha_over_k_mesh = np.linspace(alpha_over_k_min, alpha_over_k_max, num=steps)
    ax.plot(alpha_over_k_mesh, rho_func(alpha_over_k_mesh), label='$\\rho$', c=colors_theor[0])
    ax.plot(alpha_over_k_mesh, j_func(alpha_over_k_mesh), label='$j$', c=colors_theor[1])


def rho_func(aoK):
    return aoK * l / (1 + aoK * (l - 1))


def j_func(aoK):
    """
    As a function of alpha/k
    """
    up = aoK * (1 - aoK) * (1 + np.sqrt(l)) ** 2
    down = 1 + aoK * (l - 1)
    return up / down


def adjust_plot():
    plt.xlabel('Effective initiation rate $\\alpha/k$')
    plt.ylabel('$\\rho$, $j$')
    plt.ylim([-0.0, 1.02])
    plt.legend(labelspacing=0.2, frameon=False, borderaxespad=0)


def plot_j_alpha_curve(analyses=None, pdf=False):
    """
    This function plots the relation between the time-constant alpha and the normalized density rho and polymerase flows j
    """
    marker = itertools.cycle(markers)

    fig = plt.figure(num=21)
    set_figure_size(21, rows=1, page_width_frac=0.5, height_factor=height_factor)
    ax = fig.add_subplot(111)

    plot_theory(ax)

    # Add literature data
    for article in [Darzacq2007, Tantale2016]:
        # Plot density rho
        rho, rhoV, aoK = [article[key] for key in ['rho', 'rhoV', 'alpha_over_k']]
        c = 'k'
        m = next(marker)
        ax.scatter(aoK, rho, s=ms1, c=c, marker=m,  zorder=4, label=article['label'])
        ax.errorbar(aoK, rho, yerr=np.sqrt(rhoV), fmt='.', markersize=ms2, c=c,
                    elinewidth=elinewidth, alpha=alpha, zorder=3)

        # Plot flux j
        j, jV, aoK = [article[key] for key in ['j', 'jV', 'alpha_over_k']]
        ax.scatter(aoK, j, s=ms1, c=c, marker=m,  zorder=4)
        ax.errorbar(aoK, j, yerr=np.sqrt(jV), fmt='.', markersize=ms2, c=c,
                    elinewidth=elinewidth, alpha=alpha, zorder=3)
        plt.xlim([alpha_over_k_min - 1e-3, alpha_over_k_max])

    adjust_plot()
    plt.show()
    plt.tight_layout()

    # Save
    figname1 = figname + '_no_data'
    figpath = os.path.join(figures_folder, figname1)
    plt.savefig(figpath + '.png', pad_inches=0, bbox_inches='tight')
    if pdf:
        plt.savefig(figpath + '.pdf', pad_inches=0, bbox_inches='tight')

    # Add and plot our data separately for each gene
    if analyses is not None:
        num = 22
        for gene in genes:

            # Reinitialize marker and color generators
            marker = itertools.cycle(markers)
            color = itertools.cycle(colors)

            set_figure_size(num=num, rows=1, page_width_frac=0.5, height_factor=height_factor)
            fig, ax = plt.subplots(1, 1, num=num, clear=True)

            plot_theory(ax)

            # Plot different constructs separately on the same plot
            for construct_id, construct in enumerate(constructs):
                c = next(color)
                m = next(marker)
                analyses_filtered = analyses[(analyses.gene == gene)
                                             & (analyses.construct == construct)]

                # %% Load data for j
                x = aoK = analyses_filtered.loc[:, 'alpha_over_k_J']
                aoKV = analyses_filtered.loc[:, 'alpha_over_k_JV'].values
                y = j = analyses_filtered.loc[:, 'j'].values
                jV = analyses_filtered.loc[:, 'jV'].values
                xstd = np.sqrt(aoKV)
                ystd = np.sqrt(jV)

                # Plot
                ax.scatter(x, y,
                           label=f'{long_labels[construct_id]}', s=ms1, c=c, marker=m,  zorder=4)
                ax.errorbar(x, y, xerr=xstd, yerr=ystd, fmt='.', markersize=ms2, c=c,
                            elinewidth=elinewidth, alpha=alpha, zorder=3)

                # %% Load data for rho
                x = aoK = analyses_filtered.loc[:, 'alpha_over_k_rho'].values
                aoKV = analyses_filtered.loc[:, 'alpha_over_k_rhoV'].values
                y = rho = analyses_filtered.loc[:, 'rho'].values
                rhoV = analyses_filtered.loc[:, 'rhoV'].values
                xstd = np.sqrt(aoKV)
                ystd = np.sqrt(rhoV)

                # Plot
                ax.scatter(x, y,
                           s=ms1, c=c, marker=m,  zorder=4)
                ax.errorbar(x, y, xerr=xstd, yerr=ystd, fmt='.', markersize=ms2, c=c,
                            elinewidth=elinewidth, alpha=alpha, zorder=3)

            adjust_plot()
            plt.xlim([alpha_over_k_min - 1e-3, alpha_over_k_max])
            plt.title(gene_long[gene])

            plt.show()
            plt.tight_layout()

            # Save
            figname1 = figname + f'_{gene}'
            figpath = os.path.join(figures_folder, figname1)
            plt.savefig(figpath + '.png', pad_inches=0, bbox_inches='tight')
            if pdf:
                plt.savefig(figpath + '.pdf', pad_inches=0, bbox_inches='tight')
    return
