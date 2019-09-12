import itertools

import numpy as np
from matplotlib import pyplot as plt

from constants import (Darzacq2007, Tantale2016, colors_additivity, gene_long,
                       k, l)
from set_figure_size import set_figure_size


def plot_j_alpha_curve(analyses=None, pdf=False):
    """
    This function plots the relation between the time-constant and the normalized density and polymerase flows
    """

    alpha_over_k_min = 0
    alpha_over_k_max = 1 / (1 + np.sqrt(l))
    steps = 100
    ms1 = 8
    ms2 = 4
    alpha = 0.25
    height_factor = 0.65
    elinewidth = 1
    markers = ('x', 'o', '^', 'd', '*')
    # colors = ['#0072BD', '#FF7F0E']
    colors_theor = [colors_additivity['sum'],  colors_additivity['brown']]
    colors = [colors_additivity[key] for key in colors_additivity]

    marker = itertools.cycle(markers)
    long_labels = ['bac', 'no pr', 'no sh']
    gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}
    # colors = [[0, 0.251, 0], [0, 0, 1], [0.9412, 0.4706, 0], [0.502, 0.251, 0]]
    # ('b', 'brown', 'r', 'g'))

    genes = ['hb', 'sn', 'kn']
    constructs = ['bac', 'no_pr', 'no_sh']
    figname = 'j_alpha'
    # color = itertools.cycle(colors)

    fig = plt.figure(num=21)
    set_figure_size(21, rows=1, page_width_frac=0.5, height_factor=height_factor)
    ax = fig.add_subplot(111)

    def plot_theory(ax):

        alpha_over_k_mesh = np.linspace(alpha_over_k_min, alpha_over_k_max, num=steps)
        ax.plot(alpha_over_k_mesh, rho_func(alpha_over_k_mesh), label='$\\rho$', c=colors_theor[0])
        ax.plot(alpha_over_k_mesh, j_func(alpha_over_k_mesh), label='$j$', c=colors_theor[1])

    plot_theory(ax)
    # Add literature data
    for article in [Darzacq2007, Tantale2016]:
        # Plot density
        rho, rhoV, aoK = [article[key] for key in ['rho', 'rhoV', 'alpha_over_k']]
        c = 'k'  # colors[0]
        m = next(marker)
        ax.scatter(aoK, rho, s=ms1, c=c, marker=m,  zorder=4, label=article['label'])
        # xerr=np.sqrt(ktauV)
        ax.errorbar(aoK, rho, yerr=np.sqrt(rhoV), fmt='.', markersize=ms2, c=c,
                    elinewidth=elinewidth, alpha=alpha, zorder=3)

        # Plot flux
        j, jV, aoK = [article[key] for key in ['j', 'jV', 'alpha_over_k']]
        # c = next(color)
        # c = colors[1]
        # m = next(markedr)
        ax.scatter(aoK, j, s=ms1, c=c, marker=m,  zorder=4)
        # xerr=np.sqrt(ktauV)
        ax.errorbar(aoK, j, yerr=np.sqrt(jV), fmt='.', markersize=ms2, c=c,
                    elinewidth=elinewidth, alpha=alpha, zorder=3)
        plt.xlim([alpha_over_k_min - 1e-3, alpha_over_k_max])

    adjust_plot()
    plt.show()
    plt.tight_layout()
    figname1 = figname + '_no_data'
    plt.savefig(figname1 + '.png', pad_inches=0, bbox_inches='tight')
    if pdf:
        plt.savefig(figname1 + '.pdf', pad_inches=0, bbox_inches='tight')

    # Add current data
    if analyses is not None:
        # %% Plot a diagram for our data
        num = 22
        for gene_id, gene in enumerate(genes):

            marker = itertools.cycle(markers)
            color = itertools.cycle(colors)

            set_figure_size(num=num, rows=1, page_width_frac=0.5, height_factor=height_factor)
            fig, ax = plt.subplots(1, 1, num=num, clear=True)
            ax.clear()
            markersize = 10
            plot_theory(ax)

            for construct_id, construct in enumerate(constructs):

                # if gene == 'hb' and construct == 'bac':
                #     figlist = ['_without_data', '_with_data']
                # else:
                #     figlist = ['_with_data']
                # figname = '_with_data'
                # ax.clear()

                analyses_filtered = analyses[(analyses.gene == gene)
                                             & (analyses.construct == construct)]
                grouped_data = analyses_filtered.groupby(by='nc')
                # print(grouped_data)
                # if figname == '_with_data':
                # for nc in range(11, 15):

                c = next(color)
                m = next(marker)

                # j
                x = aoK = analyses_filtered.loc[:, 'alpha_over_k_J']
                # print(x)
                aoKV = analyses_filtered.loc[:, 'alpha_over_k_JV'].values
                y = j = analyses_filtered.loc[:, 'j'].values
                jV = analyses_filtered.loc[:, 'jV'].values

                # # j
                # x = aoK = grouped_data.head(1).loc[:, 'alpha_over_k_comb'].values
                # # print(x)
                # aoKV = grouped_data.head(1).loc[:, 'alpha_over_k_combV'].values
                # y = j = grouped_data.mean()['j'].values
                # jV = grouped_data.var(ddof=1).loc[:, 'j'].values
                # # print(jV)

                xstd = np.sqrt(aoKV)
                ystd = np.sqrt(jV)

                ax.scatter(x, y,
                           label=f'{long_labels[construct_id]}', s=ms1, c=c, marker=m,  zorder=4)
                ax.errorbar(x, y, xerr=xstd, yerr=ystd, fmt='.', markersize=ms2, c=c,
                            elinewidth=elinewidth, alpha=alpha, zorder=3)

                # rho
                x = aoK = analyses_filtered.loc[:, 'alpha_over_k_rho'].values
                aoKV = analyses_filtered.loc[:, 'alpha_over_k_rhoV'].values
                y = rho = analyses_filtered.loc[:, 'rho'].values
                rhoV = analyses_filtered.loc[:, 'rhoV'].values

                # # rho
                # x = aoK = grouped_data.head(1).loc[:, 'alpha_over_k_rho'].values
                # aoKV = grouped_data.head(1).loc[:, 'alpha_over_k_rhoV'].values
                # y = rho = grouped_data.mean().loc[:, 'rho'].values
                # rhoV = grouped_data.var(ddof=1).loc[:, 'rho'].values

                xstd = np.sqrt(aoKV)
                ystd = np.sqrt(rhoV)

                ax.scatter(x, y,
                           s=ms1, c=c, marker=m,  zorder=4)
                ax.errorbar(x, y, xerr=xstd, yerr=ystd, fmt='.', markersize=ms2, c=c,
                            elinewidth=elinewidth, alpha=alpha, zorder=3)

            adjust_plot()
            plt.xlim([alpha_over_k_min - 1e-3, alpha_over_k_max])
            plt.title(gene_long[gene])

            plt.show()
            plt.tight_layout()
            figname1 = figname + f'_{gene}'
            plt.savefig(figname1 + '.png', pad_inches=0, bbox_inches='tight')
            if pdf:
                plt.savefig(figname1 + '.pdf', pad_inches=0, bbox_inches='tight')

            # plot_theoretical_curve(ax, abort=False, norm=True)
            # plot_adjust(analyses_filtered)

            # ax.plot([gamma_max] * 2, ylims, '-r')
            # ax.plot([1] * 2, ylims, '-k')
            # plt.ylim(ylims)
            # plt.xlim(xlims)

            # figname = os.path.join(output_folder, 'current-density_' +
            #                        gene_long[gene] + figname)
            # fig.savefig(figname + '.png')
            # if pdf:
            #     fig.savefig(figname + '.pdf')
    return


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

    # ax.set_xscale('log')

    plt.ylim([-0.0, 1.02])
    plt.legend(labelspacing=0.2, frameon=False, borderaxespad=0)  # loc='lower left')
