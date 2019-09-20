

import matplotlib.pyplot as plt
import numpy as np

from constants import colors_additivity
from set_figure_size import set_figure_size
from theoretical_phase_parameters import (J_over_k_HD, J_over_k_LD,
                                          alpha_over_k_abortive,
                                          alpha_over_k_MC, rho_HD, rho_LD)


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
