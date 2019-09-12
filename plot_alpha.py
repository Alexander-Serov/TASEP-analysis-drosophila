
import itertools

import matplotlib.pyplot as plt
import numpy as np

from constants import alpha_additivity, colors_additivity, markers_additivity
from set_figure_size import set_figure_size
from theoretical_phase_parameters import alpha_over_k_abortive


def plot_alpha(analyses, gene='hb'):
    constructs = ['bac', 'no_pr', 'no_sh']
    # long_labels = ['no shadow', 'no primary', 'bac']
    long_labels = ['bac', 'no pr', 'no sh']
    gene_long = {'hb': 'Hunchback'}
    y_label = {'alpha_over_k_comb': 'Attempted injection rate, $\\alpha/k$'}
    x_jiggle = 0.04
    x_shifts = np.array([-1, 0, 1]) * x_jiggle
    capsize = 0

    markersize = 4
    lw = 1
    ncs = np.arange(11, 15)

    grouped_data = analyses.groupby(by=['gene', 'construct', 'nc'])
    all_means = grouped_data.mean()
    all_ps = analyses.groupby(by=['gene', 'nc']).first()
    # all_stds = grouped_data.std(ddof=1)

    quantity = 'alpha_over_k_comb'
    num = 15
    set_figure_size(num=num, rows=1, page_width_frac=0.5, clear=True, height_factor=0.7)
    fig, ax = plt.subplots(1, 1, num=num, clear=True)

    labels = [f'alpha_over_k{nc}' for nc in ncs]
    avg_data, std_data = {}, {}
    for construct in constructs:
        # filter = (analyses_norm.gene == gene) & (
        #     analyses_norm.construct == construct)
        # # plot_data = [analyses_norm[filter].loc[:, J_label].dropna().values for J_label in J_labels]
        avg_data[construct] = all_means.loc[(gene, construct, slice(None)), quantity].values
        std_data[construct] = np.sqrt(
            all_means.loc[(gene, construct, slice(None)), quantity + 'V'].values)
        # np.array(
        #     [analyses_norm[filter].loc[:, label].mean() for label in labels])
        # std_data[construct] = np.array(
        #     [np.std(analyses_norm[filter].loc[:, label], ddof=1) for label in labels])

    # TODO: remove when have actual data
    # print(avg_data)
    # avg_data['nop'] = avg_data['bac'] * 0.7
    # avg_data['nosh'] = avg_data['bac'] * 0.45
    # std_data['nop'] = std_data['nosh'] = std_data['bac'] * 0.5

    #
    avg_data['sum'] = avg_data['no_pr'] + avg_data['no_sh']
    std_data['sum'] = np.sqrt(std_data['no_pr']**2 + std_data['no_sh']**2)
    print(avg_data)
    print(std_data)
    marker_gen = itertools.cycle(markers_additivity)
    for i, construct in enumerate(constructs):
        m = next(marker_gen)
        plt.errorbar(ncs + x_shifts[i], avg_data[construct], yerr=std_data[construct],
                     fmt='-' + m,  color=colors_additivity[construct], capsize=capsize, label=long_labels[i], markersize=markersize, lw=lw)
        # plt.plot(ncs + x_shifts[i], avg_data[construct], '-o',
        #          color=colors_additivity[construct], label=long_labels[i])

    # Make a filled region for the sum
    x_fill = ncs.copy().astype(np.float)
    # print(x_fill)
    x_fill[[0, 1]] -= x_jiggle
    x_fill[[-2, -1]] += x_jiggle
    # plt.fill_between(x_fill, avg_data['sum'] - std_data['sum'], avg_data['sum'] +
    #                  std_data['sum'], facecolor=colors_additivity['sum'], alpha=alpha_additivity, edgecolor=None)

    # Add abortive theory
    abrt_limit = alpha_over_k_abortive
    # ax.plot(plt.xlim(), [abrt_limit] * 2, '--', color=colors_additivity['abrt_limit'], lw=1)

    # Add p-values
    p = all_ps[quantity + '_p_value'].values
    y = np.max(avg_data[construct] + std_data[construct]) * 1.1
    for i, nc in enumerate(ncs):
        x = nc

        if not np.isnan(p[i]):
            ax.text(x, y, f'p={p[i]:.2f}', rotation=90, va='bottom', ha='center')

    # ax.boxplot(plot_data, positions=ncs, showfliers=True)
    # Mean line
    # ax.plot(ncs, avg_data)
    plt.ylim(ymin=0)
    # _, ymax = plt.ylim()
    # print(y, y * 1.5)
    plt.ylim(ymax=y * 1.4)

    plt.xlabel('Nuclear cycle')
    plt.ylabel(y_label[quantity])
    plt.xticks(ncs)
    plt.title(gene_long[gene])
    # plt.legend(loc='upper left')
    plt.tight_layout()
    plt.show()

    figname = 'additivity_alpha_' + gene
    fig.savefig(figname + '.pdf', pad_inches=0, bbox_inches='tight')
    fig.savefig(figname + '.png', pad_inches=0, bbox_inches='tight')
