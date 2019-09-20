"""
Perform Welch's unequal variances t-test to evaluate whether alpha is the same across genes, constructs or ncs.
Each of the three parts of the code corresponds to one changing parameter.

All of the tests are performed for alpha since we saw that it is close to being normally distributed, while tau isn't.
Similarities between other parameters are tested for historical reasons and are not re-used later on.
"""

import numpy as np
import pandas as pd
from scipy.stats import t as student

from constants import gene_labels


def perform_welchs_test(analyses_in):
    analyses = analyses_in.copy()

    genes = set(analyses.gene)
    constructs = set(analyses.construct)
    gene_ids = set(analyses.gene_id)
    ncs = set(analyses.index.get_level_values(1))
    quantities = ['r', 'j', 'alpha_over_k_comb', 'tau', 'alpha_comb']

    grouped = analyses.groupby(by=['gene', 'nc', 'construct'])
    idx = pd.IndexSlice

    # Initialize a temporary data frame
    calc = pd.DataFrame()
    for col in ['t', 'nu']:
        calc[col] = np.nan

    # %% Test for similarities across constructs
    # Some tests are only performed if all three constructs are availables
    has_all_constructs = 'bac' in constructs and 'no_pr' in constructs and 'no_sh' in constructs
    if not has_all_constructs:
        print(">> Skipping Welch's test across constructs: some of the constructs were not detected <<")
    else:
        print('>> p-values across constructs <<')
        for quantity in quantities:
            # Load/calculate means, variances and number of samples
            if quantity in ['r', 'j']:
                calc['means'] = grouped[quantity].mean()
                calc['vars'] = grouped[quantity].apply(
                    lambda group: np.var(group[~np.isnan(group)], ddof=1))

            elif quantity in ['alpha_over_k_comb', 'tau', 'alpha_comb']:
                calc['means'] = grouped[quantity].first()
                calc['vars'] = grouped[quantity + 'V'].first()

            calc['n'] = grouped[quantity].count()

            # Calculate independently for each gene and nc
            for gene in genes:
                for nc in ncs:
                    bac = calc.loc[(gene, nc, 'bac'), :]
                    no_sh = calc.loc[(gene, nc, 'no_sh'), :]

                    p, t, nu = welchs_test(
                        x1Mean=no_sh.means, x1V=no_sh.vars,
                        n1=no_sh.n, x2Mean=bac.means, x2V=bac.vars, n2=bac.n)

                    calc.loc[(gene, nc, 'bac'), ['t']] = t
                    calc.loc[(gene, nc, 'bac'), ['nu']] = nu
                    label = quantity + '_p_value'
                    calc.loc[(gene, nc, 'bac'), label] = p

                    # Copy the p-value back into the analyses table
                    gene_filter = analyses.gene == gene
                    analyses.loc[idx[gene_filter, nc], label] = p

                    # Print out results for `alpha_comb`
                    if not np.isnan(t) and quantity is 'alpha_comb':
                        print(f'{gene}, nc{nc}, bac, no_sh:\tp for alpha_comb: {p:.2g}')

    # %% Test across ncs for only alpha similarities.
    # Basically I compare nc14 to all other ncs for each gene-construct combination
    print('\n>> p-values across ncss <<')
    quantity = 'alpha_comb'
    nc1 = 14

    # Load means, vars and sample numbers
    calc['means'] = grouped[quantity].first()
    calc['vars'] = grouped[quantity + 'V'].first()
    calc['n'] = grouped[quantity].count()

    for gene in genes:
        for construct in constructs:
            for nc2 in ncs:
                if nc2 >= nc1:
                    continue
                data1 = calc.loc[(gene, nc1, construct), :]
                data2 = calc.loc[(gene, nc2, construct), :]

                p, t, nu = welchs_test(
                    x1Mean=data1.means, x1V=data1.vars,
                    n1=data1.n, x2Mean=data2.means, x2V=data2.vars, n2=data2.n)

                print(f'{gene}, {construct}, nc{nc2}, nc{nc1}:\tp for alpha: {p:.2g}')

    # %% Similarity test for alpha across genes
    if len(genes) < 2:
        print("\n>> Skipping Welch's test across genes: not enough genes <<")
    else:
        print('\n>> p-values across genes <<')

        grouped = analyses.groupby(by=['gene', 'nc', 'construct'])
        index = pd.MultiIndex.from_product([ncs, genes, genes])

        quantities = ['alpha_comb', 'tau']
        print('p-values across genes: ')
        for construct in constructs:
            for nc in ncs:
                for gene_id1 in gene_ids:
                    for gene_id2 in gene_ids:
                        if gene_id2 <= gene_id1:
                            continue
                        p = {}
                        for quantity in quantities:

                            # Create a temporary data frame
                            calc = pd.DataFrame()
                            calc['means'] = grouped[quantity].first()
                            calc['vars'] = grouped[quantity + 'V'].first()
                            calc['n'] = grouped[quantity].count()

                            gene1, gene2 = [gene_labels[i] for i in [gene_id1, gene_id2]]
                            data1 = calc.loc[(gene1, nc, construct), :]
                            data2 = calc.loc[(gene2, nc, construct), :]

                            p[quantity], t, nu = welchs_test(
                                x1Mean=data1.means, x1V=data1.vars,
                                n1=data1.n, x2Mean=data2.means, x2V=data2.vars, n2=data2.n)
                        print(
                            f'{construct}, nc{nc}, {gene1}, {gene2}:\tp for {quantities[0]} - {p[quantities[0]]:.2g},\tp for {quantities[1]} - {p[quantities[1]]:.2g}')
    return analyses


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


def welchs_test_arrays(x1, x2):
    """
    A spearate implementation that automatically calculates means, variances and the number of samples from the input lists of values.

    Not currently used in the code because the number of samples is calculated earlier, when the estimators were combined.
    It is kept as a separate number for the mixed estimators.
    """

    x1, x2 = x1[~np.isnan(x1)], x2[~np.isnan(x2)]
    n1, n2 = len(x1), len(x2)

    x1Mean = np.mean(x1)
    x2Mean = np.mean(x2)

    x1V = np.var(x1, ddof=1)
    x2V = np.var(x2, ddof=1)

    return welchs_test(x1Mean, x1V, n1, x2Mean, x2V, n2)
