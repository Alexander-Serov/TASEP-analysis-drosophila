
import itertools
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constants import LV, L, colors_additivity, k, kV, l, lV
from least_squares_BCES import least_squares_BCES
from plot_theory_current_density import plot_theoretical_curve
from reinit_folder import reinit_folder
from set_figure_size import set_figure_size

height_factor = 0.65

output_folder = 'current-density'


def plot_current_density_diagram(analyses_in, I_est, num, I_estV=0, pdf=True):

    analyses = analyses_in.copy()
    long_labels = ['bac', 'no pr', 'no sh']
    gene_long = {'hb': 'Hunchback'}
    # colors = [[0, 0.251, 0], [0, 0, 1], [0.9412, 0.4706, 0], [0.502, 0.251, 0]]
    # ('b', 'brown', 'r', 'g'))

    genes = ['hb']
    constructs = ['bac', 'no_pr', 'no_sh']
    markers = ('x', 'o', '^', 'd', '*')
    # colors = ('#F16822', '#68972F', '#BE2026', 'mediumseagreen')
    colors = [colors_additivity[key] for key in colors_additivity]
    ms1 = 8
    ms2 = 3
    alpha = 0.25
    elinewidth = 1

    reinit_folder(output_folder)

    iterables = pd.MultiIndex.from_product([genes, constructs])
    Ts = pd.DataFrame(columns=['T', 'TV', 'Tstd', 'L', 'LV', 'Lstd'], index=iterables)

    # %% Plot a diagram with literature data
    fig = set_figure_size(num=num, rows=1, page_width_frac=0.5, height_factor=height_factor)
    ax = fig.subplots(1, 1)
    ax.clear()
    marker = itertools.cycle(markers)
    color = itertools.cycle(colors)
    Darzacq2007 = {'rho': 0.023,
                   'rhostd': 0.008,
                   'JoK': 1.3e-5,
                   'JoKstd': 0.4e-5,
                   'label': 'Darzacq2007'}

    Tantale2016 = {'rho': 0.19,
                   'rhostd': 0.17,
                   'JoK': 3.9e-3,
                   'JoKstd': 3.2e-3,
                   'label': 'Tantale2016'}

    # Plot
    for data in [Darzacq2007, Tantale2016]:
        c = next(color)
        m = next(marker)
        x, xstd, y, ystd = [data[key] for key in ['rho', 'rhostd', 'JoK', 'JoKstd']]
        ax.scatter(x, y, s=ms1, c=c, marker=m,  zorder=4, label=data['label'])
        ax.errorbar(x, y, xerr=xstd, yerr=ystd, fmt='.', markersize=ms2, c=c,
                    elinewidth=elinewidth, alpha=alpha, zorder=3)
    # plt.ylim(ymin=-1e-3)
    plot_theoretical_curve(ax, abort=False)

    # Adjust & save
    def plot_adjust(analyses_filtered=None):
        plt.xlabel('Site occupation density, $\\rho$')
        plt.ylabel('Polymerase flux, $J/k$')
        plt.xlim([0, 1.02])
        plt.ylim(-2e-4)
        plt.yticks(np.arange(0, 0.018, step=0.005))

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
    # plt.ylim(ymin=-1e-3)

    figname = os.path.join(output_folder, 'current-density_literature')
    fig.savefig(figname + '.png')
    if pdf:
        fig.savefig(figname + '.pdf')

    # %% Plot a diagram for our data
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

            x = rho = analyses_filtered.loc[:, 'rho'].values
            rhoV = analyses_filtered.loc[:, 'rhoV'].values
            y = JoK = analyses_filtered.loc[:, 'JoK'].values
            JoKV = analyses_filtered.loc[:, 'JoKV'].values

            xstd = np.sqrt(rhoV)
            ystd = np.sqrt(JoKV)

            c = next(color)
            m = next(marker)
            ax.scatter(x, y,
                       label=f'{long_labels[construct_id]}', s=ms1, c=c, marker=m,  zorder=4)
            ax.errorbar(x, y, xerr=xstd, yerr=ystd, fmt='.', markersize=ms2, c=c,
                        elinewidth=elinewidth, alpha=alpha, zorder=3)
            # else:
            #     # For the plot without data, plot data from other articles, format: (rho, J/k)
            #     nMRNA = np.array([200, 400])
            #     nGenes = 200
            #     Darzacq2007 = {'rho': (nMRNA * 0.7 * l / 2300 / nGenes).mean(),
            #                    'rhostd': np.std(nMRNA * 0.7 * l / 2300 / nGenes, ddof=1),
            #                    'JoK': (nMRNA * 0.39 * 0.00159 / nGenes / 2300 * 32).mean(),
            #                    'JoKstd': np.std(nMRNA * 0.39 * 0.00159 / nGenes / 2300 * 32, ddof=1)}
            #
            #     kElong = 4e3 / 60  # bp/s
            #     Tantale2016 = {'rho': l / (4.1 * kElong),
            #                    'rhostd': 0,
            #                    'JoK': 1 / 4.1 / kElong,
            #                    'JoKstd': 0}
            #     ax.errorbar(Darzacq2007['rho'], Darzacq2007['JoK'], xerr=Darzacq2007['rhostd'], yerr=Darzacq2007['JoKstd'], fmt='x', c='r',  # markersize=ms2, c=c,
            #                 elinewidth=1, alpha=1, zorder=3)  # , label='Darzacq2007')
            #
            #     ax.errorbar(Tantale2016['rho'], Tantale2016['JoK'], xerr=Tantale2016['rhostd'], yerr=Tantale2016['JoKstd'], fmt='x', c='g',  # markersize=ms2, c=c,
            #                 elinewidth=1, alpha=1, zorder=3)  # , label='Tantale2016')

            plot_theoretical_curve(ax, abort=False)
            plot_adjust(analyses_filtered)

            # ax.plot([gamma_max] * 2, ylims, '-r')
            # ax.plot([1] * 2, ylims, '-k')

            figname = os.path.join(output_folder, 'current-density_' +
                                   gene_long[gene] + figname)
            fig.savefig(figname + '.png', pad_inches=0, bbox_inches='tight')
            if pdf:
                fig.savefig(figname + '.pdf', pad_inches=0, bbox_inches='tight')

    return analyses
