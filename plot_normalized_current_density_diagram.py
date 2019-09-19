"""
Plot the normalized current-density diagram with and without earlier literature data.
"""


import itertools
import os

import matplotlib.pyplot as plt
import numpy as np

from constants import (LV, Darzacq2007, L, Tantale2016, colors_additivity,
                       figures_folder)
from plot_theory_current_density import plot_theoretical_curve
from set_figure_size import set_figure_size
from theoretical_phase_parameters import (J_over_k_HD, J_over_k_LD,
                                          alpha_over_k_MC, rho_HD, rho_LD)

height_factor = 0.65
ylims = [0, 1.01]
xlims = [0, 1.02]
gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}


def plot_normalized_current_density_diagram(analyses_in, num, pdf=True):

    # Constants
    markers = ('x', 'o', '^', 'd', '*')
    colors = [colors_additivity[key] for key in colors_additivity]
    ms1 = 8  # marker size
    ms2 = 3  # smaller markers for the error bars
    alpha = 0.25
    elinewidth = 1  # line width for the errorbars
    long_labels = ['bac', 'no pr', 'no sh']

    analyses = analyses_in.copy()
    genes = set(analyses.gene)
    constructs = set(analyses.construct)

    # %% Plot a diagram with literature data
    fig = set_figure_size(num=num, rows=1, page_width_frac=0.5,
                          height_factor=height_factor, clear=True)
    ax = fig.subplots(1, 1)
    marker = itertools.cycle(markers)
    color = itertools.cycle(colors)

    # Plot
    for data in [Darzacq2007, Tantale2016]:
        c = 'k'  # color
        m = next(marker)
        rho, rhoV, j, jV = [data[key] for key in ['rho', 'rhoV', 'j', 'jV']]
        ax.scatter(rho, j, s=ms1, c=c, marker=m,  zorder=4, label=data['label'])
        ax.errorbar(rho, j, xerr=rhoV**(1 / 2), yerr=jV**(1 / 2), fmt='.', markersize=ms2, c=c,
                    elinewidth=elinewidth, alpha=alpha, zorder=3)

    # Add the current-density curve
    plot_theoretical_curve(ax)
    plot_adjust()
    plt.ylim(ymin=-0.02)

    # Save figure
    figname = os.path.join(figures_folder, 'current-density_literature')
    fig.savefig(figname + '.png', pad_inches=0, bbox_inches='tight')
    if pdf:
        fig.savefig(figname + '.pdf', pad_inches=0, bbox_inches='tight')

    # %% Plot a diagram with our data
    num += 1
    for gene in genes:
        marker = itertools.cycle(markers)
        color = itertools.cycle(colors)

        set_figure_size(num=num, rows=1, page_width_frac=0.5, height_factor=height_factor)
        fig, ax = plt.subplots(1, 1, num=num, clear=True)

        # Plot data grouped by gene and construct
        for construct_id, construct in enumerate(constructs):
            figname = '_with_data'
            c = next(color)
            m = next(marker)

            analyses_filtered = analyses[
                (analyses.gene == gene) & (analyses.construct == construct)]

            x = analyses_filtered.loc[:, 'rho'].values
            rV = analyses_filtered.loc[:, 'rhoV'].values
            y = analyses_filtered.loc[:, 'j'].values
            jV = analyses_filtered.loc[:, 'jV'].values

            xstd = np.sqrt(rV)
            ystd = np.sqrt(jV)

            ax.scatter(x, y,
                       label=f'{long_labels[construct_id]}', s=ms1, c=c, marker=m,  zorder=4)
            ax.errorbar(x, y, xerr=xstd, yerr=ystd, fmt='.', markersize=ms2, c=c,
                        elinewidth=elinewidth, alpha=alpha, zorder=3)

            # Add the current-density diagram
            plot_theoretical_curve(ax)

            # Adjust the plot
            plot_adjust()
            plt.title(f'{gene_long[gene]}')

            # Add fits of the current-density diagram: transit time and the effective gene length
            T, TV, L, LV = analyses_filtered[['T', 'TV', 'L', 'LV']].values[0]
            L_rnd, Lstd_rnd = np.round(np.array([L, np.sqrt(LV)]) / 1e3, decimals=2)
            ax.text(
                0.57, 0.05, f'$T = {T:.2f} \pm {np.sqrt(TV):.2f}$ min\n$L = {L_rnd:.1f} \pm {Lstd_rnd:.1f}$ kb', transform=ax.transAxes, va='bottom', ha='left', fontsize=8, weight='normal', style='normal', family='sans-serif')

            # Save figure
            figname = os.path.join(figures_folder, 'current-density_' +
                                   gene_long[gene] + figname)
            fig.savefig(figname + '.png', pad_inches=0, bbox_inches='tight')
            if pdf:
                fig.savefig(figname + '.pdf', pad_inches=0, bbox_inches='tight')

    return


def plot_adjust():
    """
    Adjust the plot
    """
    plt.xlabel('Site occupation density $\\rho$')
    plt.ylabel('Normalized flux $j$')
    plt.xlim(xlims)
    plt.ylim(ylims)

    plt.legend(loc='upper left', labelspacing=0.2, frameon=False, borderaxespad=0)
    plt.tight_layout()


def plot_theoretical_curve(ax):
    """
    Plot current vs density for the whole available range of densities.
    This inlcudes LD, MC and HD phases.
    Formulas based on Shaw2003 and the accompanying manuscript.

    Notation:

    rho                 # site occupation density
    rho/l               # polymerase density
    J                   # polymerase current
    alphaT = alpha/k    # Dimensionless injection attempt rate \tilde{alpha}
    betaT = beta/k    # Dimensionless exit attempt rate \tilde{beta}
    J = sT = s/k    # particle current
    alphaTMC = betaTMC # Dimensionless rate for the max. current regime
    """

    steps = 100
    colorLD = colors_additivity['sum']  # '#0072BD'
    colorHD = colors_additivity['brown']  # '#FF7F0E'

    a_mesh = np.linspace(0, alpha_over_k_MC, num=steps, endpoint=True)
    b_mesh = a_mesh

    JoK_MC = J_over_k_LD(alpha_over_k_MC)
    J_fact = JoK_MC

    ax.plot(rho_LD(a_mesh), J_over_k_LD(a_mesh) / J_fact, c=colorLD)
    ax.plot(rho_HD(b_mesh), J_over_k_HD(b_mesh) / J_fact, c=colorHD)

    return
