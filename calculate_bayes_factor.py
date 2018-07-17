

import numpy as np
from numpy import log
from numpy import pi
from scipy.special import gammaln


def calculate_bayes_factor(dataset1, dataset2, n_pi, mu_pi, V_pi):
    """
    Calculate the Bayes factor for the evidence that two datasets should be described by two different Gaussian distributions, as compared to a single Gaussian.
    Basically it answers the question of whether the two distributions are different.
    """

    def pre_process(dataset):
        # dataset = np.asarray(dataset)
        dataset = dataset[~np.isnan(dataset)]
        return dataset

    dataset1 = pre_process(dataset1)
    dataset2 = pre_process(dataset2)

    dataset = np.concatenate([dataset1, dataset2])
    n1 = len(dataset1)
    n2 = len(dataset2)
    n = n1 + n2
    # print([n1, n2, n])

    # Calculate mean and biased variances
    s1 = dataset1.mean()
    s2 = dataset2.mean()
    s = dataset.mean()
    V1 = dataset1.var()
    V2 = dataset2.var()
    V = dataset.var()
    std1 = dataset1.std()
    std2 = dataset2.std()
    # print([s1, s2, s, V1, V2, V, n_pi, mu_pi, V_pi])

    log_A0 = (log(2) + 1 / 2 * n_pi + (n_pi - 1) / 2 * log(pi) +
              (n_pi - 3) / 2 * log(n_pi * V_pi) - gammaln((n_pi - 3) / 2))
    # print(log_A0)

    def log_prob(s, n, V):
        return (
            log_A0
            + (1 - n - n_pi) / 2 * log(pi)
            - log(2)
            - 1 / 2 * log(n + n_pi)
            + gammaln((n + n_pi - 3) / 2)
            - (n + n_pi - 3) / 2
            * log(n * V + n_pi * V_pi + n * n_pi / (n + n_pi) * (s - mu_pi) ** 2)
        )

    log_B = log_prob(s1, n1, V1) + log_prob(s2, n2, V2) - log_prob(s, n, V)

    log10_B = log_B / log(10)
    means = [s1, s2]
    stds = [std1, std2]

    return [log10_B, means, stds]
