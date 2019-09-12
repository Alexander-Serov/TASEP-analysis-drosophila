# plot_theory_current_density
"""
Plot current vs density for the whole range of densities. This will inlcude LD, MC and HD phases.
Formulas based on Shaw2003

rho                 # site occupation density
rho/l               # polymerase density
J                   # polymerase current
alphaT = alpha/k    # Dimensionless injection attempt rate \tilde{alpha}
betaT = beta/k    # Dimensionless exit attempt rate \tilde{beta}
# Dimensionless polymerase (particle, not bp) injection rate. This is just dimensionless current J
J = sT = s/k    # particle current
alphaTMC = betaTMC # Dimensionless rate for the max. current regime

As long as I'm in LD or HD, the current and density both depend on alpha (or beta) only. So I can just independently plot two parts of the parametric curve
"""

# if __name__ == '__main__':
#     try:
#         has_run
#     except NameError:
#         %matplotlib tk
#         %load_ext autoreload
#         %autoreload 2
#         has_run = 1
#     else:
#         print("Graphic interface NOT re-initialized")

import matplotlib.pyplot as plt
import numpy as np

from constants import colors_additivity, k, l, tau_rev
from set_figure_size import set_figure_size
from theoretical_phase_parameters import (J_over_k_HD, J_over_k_LD,
                                          alpha_over_k_abortive,
                                          alpha_over_k_MC, rho_HD, rho_LD)


def plot_theoretical_curve(ax, abort=True, norm=False):
    """
    abort --- whether the abortive theory prediction should be plotted
    """

    # Plot

    steps = 100
    da = alpha_over_k_MC / steps
    a_mesh = np.linspace(0, alpha_over_k_MC, num=steps, endpoint=True)
    b_mesh = a_mesh
    colorLD = colors_additivity['sum']  # '#0072BD'
    colorHD = colors_additivity['brown']  # '#FF7F0E'

    rho_MC = rho_LD(alpha_over_k_MC)
    JoK_MC = J_over_k_LD(alpha_over_k_MC)

    if norm:
        # rho_fact = rho_MC
        J_fact = JoK_MC
    else:
        J_fact = 1
    # print('Hello!', norm)

    ax.plot(rho_LD(a_mesh), J_over_k_LD(a_mesh) / J_fact, c=colorLD)
    ax.plot(rho_HD(b_mesh), J_over_k_HD(b_mesh) / J_fact, c=colorHD)

    # Add abortive values
    if abort:
        lw = 0.75
        abrt_color = 'k'
        rho_abrt = rho_LD(alpha_over_k_abortive)
        J_abrt = J_over_k_LD(alpha_over_k_abortive)
        ymin, _ = plt.ylim()
        ax.plot([0, rho_abrt], [J_abrt] * 2, '--', lw=lw, color=abrt_color)
        ax.plot([rho_abrt] * 2, [ymin, J_abrt], '--', lw=lw, color=abrt_color)

    return


if __name__ == '__main__':
    height_factor = 3
    num = 11
    set_figure_size(num=num, rows=1, page_width_frac=0.5, height_factor=height_factor)
    fig, ax = plt.subplots(1, 1, num=num)
    plot_theoretical_curve(ax)

    plt.xlabel('rho')
    plt.ylabel('J')
    plt.xlim(0)
    plt.ylim(0)

    plt.show()
