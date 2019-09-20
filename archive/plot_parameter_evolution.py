

import itertools
import os

import matplotlib.pyplot as plt
import numpy as np

from constants import colors_additivity, figures_folder, markers_additivity
from set_figure_size import set_figure_size


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
            set_figure_size(num=num, rows=1, page_width_frac=0.5, clear=True, height_factor=0.7)
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
