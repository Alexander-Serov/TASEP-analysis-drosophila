"""
Perform Welch's unequal variances t-test to evaluate whenther the additivity hypothesis can be refuted for the given data sets.
The code calculates the characteristics of the random variable equal to the sum of values of the no primary and no shadow constructs and compares it to the distribution of bac.

Another part of the code checks whether the bac value is the same across different constructs.
"""

import numpy as np
import pandas as pd
from scipy.stats import t as student

from constants import tau_rev, tau_rev_n, tau_revV


def perform_additivity_test(analyses_in):

    analyses = analyses_in.copy()
    quantities = ['r', 'j', 'alpha_over_k_comb', 'tau', 'alpha_comb']
    genes = set(analyses.gene)
    ncs = set(analyses.index.get_level_values(1))
    idx = pd.IndexSlice

    grouped = analyses.groupby(by=['gene', 'nc', 'construct'])

    # Evaluate additivity
    for quantity in quantities:
        calc = pd.DataFrame()

        # A special condition for alpha since we don't have the distribution
        if quantity in ['r', 'j']:
            calc['means'] = grouped[quantity].mean()
            calc['vars'] = grouped[quantity].apply(
                lambda group: np.var(group[~np.isnan(group)], ddof=1))

        elif quantity in ['alpha_over_k_comb', 'tau', 'alpha_comb']:
            calc['means'] = grouped[quantity].first()
            calc['vars'] = grouped[quantity + 'V'].first()

        calc['n'] = grouped[quantity].count()

        # Add sum of constructs
        for col in ['t', 'nu']:
            calc[col] = np.nan

        for gene in genes:
            for nc in ncs:
                bac = calc.loc[(gene, nc, 'bac'), :]
                no_pr = calc.loc[(gene, nc, 'no_pr'), :]
                no_sh = calc.loc[(gene, nc, 'no_sh'), :]

                p, t, nu = welchs_test(x1Mean=no_sh.means, x1V=no_sh.vars,
                                       n1=no_sh.n, x2Mean=bac.means, x2V=bac.vars, n2=bac.n)

                calc.loc[(gene, nc, 'bac'), ['t']] = t
                calc.loc[(gene, nc, 'bac'), ['nu']] = nu
                label = quantity + '_p_value'
                calc.loc[(gene, nc, 'bac'), label] = p

                # Copy the p-value back into the analyses table
                gene_filter = analyses.gene == gene
                analyses.loc[idx[gene_filter, nc], label] = p

                # p - value for additivity
                if not np.isnan(t):
                    p = 2 * (student.cdf(-np.abs(t), df=nu))
                    label = quantity + '_p_value'
                    # Copy the p-value back into the analyses table
                    calc.loc[(gene, nc, 'bac'), label] = p

                    # Copy the p-value back into the analyses table
                    gene_filter = analyses.gene == gene
                    analyses.loc[idx[gene_filter, nc], label] = p

    return analyses


def welchs_test_arrays(x1, x2):
    """
    A spearate implementation not currently used in the code
    """

    x1, x2 = x1[~np.isnan(x1)], x2[~np.isnan(x2)]
    n1, n2 = len(x1), len(x2)

    x1Mean = np.mean(x1)
    x2Mean = np.mean(x2)

    x1V = np.var(x1, ddof=1)
    x2V = np.var(x2, ddof=1)

    return welchs_test(x1Mean, x1V, n1, x2Mean, x2V, n2)


def welchs_test(x1Mean, x1V, n1, x2Mean, x2V, n2):
    """
    Can accept arrays as input
    """
    # Test statistics
    t = (x1Mean - x2Mean) / \
        np.sqrt(x1V / n1 + x2V / n2)
    nu = ((x1V / n1 + x2V / n2)**2
          / (x1V**2 / n1 / (n1 - 1) + x2V**2 / n2 / (n2 - 1))
          )
    # p-value for additivity
    if not np.isnan(t):
        p = 2 * (student.cdf(-np.abs(t), df=nu))
    else:
        p = np.nan

    return p, t, nu
