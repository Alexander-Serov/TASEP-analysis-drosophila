import itertools
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt

from constants import (LV, Darzacq2007, L, Tantale2016, colors_additivity,
                       figures_folder, gene_long, l, markers_additivity)
from support import (J_over_k_HD, J_over_k_LD, alpha_over_k_abortive,
                     alpha_over_k_MC, rho_HD, rho_LD, set_figure_size)

# Constants
alpha_over_k_min = 0
alpha_over_k_max = 1 / (1 + np.sqrt(l))  # theoretical maximum
colors_theor = [colors_additivity['sum'],  colors_additivity['brown']]


def plot_theory(ax):
    """Plot theoretical curves for rho and j"""
    steps = 100

    alpha_over_k_mesh = np.linspace(alpha_over_k_min, alpha_over_k_max, num=steps)
    ax.plot(alpha_over_k_mesh, rho_func(alpha_over_k_mesh), label='$\\rho$', c=colors_theor[0])
    ax.plot(alpha_over_k_mesh, j_func(alpha_over_k_mesh), label='$j$', c=colors_theor[1])

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
    plt.ylim([-0.0, 1.02])
    plt.legend(labelspacing=0.2, frameon=False, borderaxespad=0)


def plot_j_alpha_curve(analyses=None, pdf=False):
    """
    This function plots the relation between the time-constant alpha and the normalized density rho and polymerase flows j
    """
    # Constants

    ms1 = 8
    ms2 = 4  # size of the marker for errorbar plot to superimpose on the the scatter data plot
    alpha = 0.25
    height_factor = 0.65
    elinewidth = 1  # errorbar line width
    markers = ('x', 'o', '^', 'd', '*')

    colors = [colors_additivity[key] for key in colors_additivity]
    long_labels = ['bac', 'no pr', 'no sh']
    gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}
    genes = ['hb', 'sn', 'kn']
    constructs = ['bac', 'no_pr', 'no_sh']
    figname = 'j_alpha'

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


def plot_normalized_current_density_diagram(analyses_in, num, pdf=True):
    """
    Plot the normalized current-density diagram with and without earlier literature data.
    """

    # Constants
    height_factor = 0.65
    ylims = [0, 1.01]
    xlims = [0, 1.02]
    gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}

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
    plt.xlim(xlims)
    plt.ylim(ylims)
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
            plt.xlim(xlims)
            plt.ylim(ylims)
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

    plt.legend(loc='upper left', labelspacing=0.2, frameon=False, borderaxespad=0)
    plt.tight_layout()

    return


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


def plot_parameter_evolution(analyses, pdf=False):
    """
    Plot changes in j, rho, alpha and tau across ncs for different genes and constructs.
    """
    ncs = np.arange(11, 15)
    genes = set(analyses.gene)
    constructs = set(analyses.construct)
    long_labels = {'bac': 'bac', 'no_pr': 'no pr', 'no_sh': 'no sh'}
    gene_long = {'hb': 'hunchback', 'kn': 'knirps', 'sn': 'snail'}
    y_label = {'j': 'Normalized flux $j$',
               'rho': 'Site occupation density $\\rho$', 'tau': 'Residence time $\\tau$ (s)', 'alpha_comb': 'Initiation rate $\\alpha$ (pol/min)'}

    # Add extra jiggle to be able to distinguish overlapping data points
    x_jiggle = 0.04
    x_shifts = np.array([-1, 0, 1]) * x_jiggle

    # Plot parameters
    capsize = 0
    markersize = 4
    lw = 1  # line width

    for gene in genes:
        grouped_data = analyses.groupby(by=['gene', 'construct', 'nc'])
        all_means = grouped_data.mean()
        all_stds = grouped_data.std(ddof=1)
        all_ps = analyses.groupby(by=['gene', 'nc']).first()

        for quantity in ['j', 'rho', 'tau', 'alpha_comb']:
            ymaxs = {'j': 0.36, 'rho': 0.27, 'tau': 103, 'alpha_comb': 12}
            num = 12
            set_figure_size(num=num, rows=1, page_width_frac=0.5,
                            clear=True, height_factor=0.7)
            fig, ax = plt.subplots(1, 1, num=num, clear=True)
            avg_data, std_data = {}, {}
            for construct in constructs:
                if quantity in ['rho', 'j']:
                    avg_data[construct] = all_means.loc[(
                        gene, construct, slice(None)), quantity].values
                    std_data[construct] = all_stds.loc[(
                        gene, construct, slice(None)), quantity].values

                elif quantity in ['tau', 'alpha_comb']:
                    avg_data[construct] = all_means.loc[(
                        gene, construct, slice(None)), quantity].values
                    std_data[construct] = np.sqrt(
                        all_means.loc[(gene, construct, slice(None)), quantity + 'V'].values)

            # Prepare a marker generator and plot the data with errorbars
            marker_gen = itertools.cycle(markers_additivity)
            for i, construct in enumerate(constructs):
                m = next(marker_gen)
                plt.errorbar(
                    ncs + x_shifts[i], avg_data[construct],
                    yerr=std_data[construct],
                    fmt='-' + m,  color=colors_additivity[construct],
                    capsize=capsize, label=long_labels[construct],
                    markersize=markersize, lw=lw)

            # Adjust plot
            plt.xlabel('Nuclear cycle')
            plt.ylabel(y_label[quantity])
            plt.ylim(ymin=0, ymax=ymaxs[quantity])

            plt.xticks(ncs)
            plt.title(gene_long[gene])

            plt.tight_layout()
            plt.show()

            # Save figure
            figname = 'additivity_' + quantity + '_' + gene
            figpath = os.path.join(figures_folder, figname)
            fig.savefig(figpath + '.png', pad_inches=0, bbox_inches='tight')
            if pdf:
                fig.savefig(figpath + '.pdf', pad_inches=0, bbox_inches='tight')

    return


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

    a_mesh = np.linspace(0, alpha_over_k_MC, num=steps, endpoint=True)
    b_mesh = a_mesh
    colorLD = colors_additivity['sum']  # '#0072BD'
    colorHD = colors_additivity['brown']  # '#FF7F0E'

    JoK_MC = J_over_k_LD(alpha_over_k_MC)

    J_fact = JoK_MC

    ax.plot(rho_LD(a_mesh), J_over_k_LD(a_mesh) / J_fact, c=colorLD)
    ax.plot(rho_HD(b_mesh), J_over_k_HD(b_mesh) / J_fact, c=colorHD)

    return


if __name__ == '__main__':
    """Some testing"""
    height_factor = 3
    num = 11
    set_figure_size(num=num, rows=1, page_width_frac=0.5, height_factor=height_factor)
    _, ax = plt.subplots(1, 1, num=num)
    plot_theoretical_curve(ax)

    plt.xlabel('rho')
    plt.ylabel('J')
    plt.xlim(0)
    plt.ylim(0)

    plt.show()
