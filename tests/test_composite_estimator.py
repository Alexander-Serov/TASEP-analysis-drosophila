"""

isort:skip_file
"""
#
# try:
#     has_run
# except NameError:
#     %matplotlib tk
#     %load_ext autoreload
#     %autoreload 2
#     has_run = 1
# else:
#     print("Graphic interface NOT re-initialized")

import numpy as np
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from BLUE_estimator import BLUE_estimator


def test_composite_estimator():

    x1 = 87.8
    x2 = 89.5
    data = np.array([[88.0, 87.5],
                     [82.5, 80.0],
                     [83.0, 86.5],
                     [73.5, 79.0],
                     [79, 84.5],
                     [82, 83.5],
                     [83, 79.8],
                     [80.8, 84],
                     [81, 83],
                     [79, 79],
                     [64, 76],
                     [80.5, 83.8],
                     [83, 87],
                     [81.5, 78.5]
                     ])

    experts = np.array([89.5, 82.5, 85.8, 76.3, 83.3, 83.8,
                        85.0, 81.3, 81.8, 81.0, 67.5, 83, 85, 82])
    true_yield = np.array([87.8, 87.3, 85.3, 76.8, 78.3, 89.0,
                           82.5, 84, 82.3, 80.8, 68.3, 83, 85, 81.8])
    np.sqrt(experts.var())
    np.std(true_yield - experts, ddof=0)
    np.std(true_yield - data[:, 0], ddof=0)
    np.std(true_yield - data[:, 1], ddof=0)
    # x1V = np.var(data[:, 0], ddof=1)  # 3.92**2
    # x2V = np.var(data[:, 1], ddof=1)  # 3.06**2
    # ddof is 0 because it's a maximum likelihood estimator for the normal (Wishart) distribution
    cov_matrix = np.cov(data, rowvar=False, ddof=0)
    np.sqrt(cov_matrix)

    d1 = data.T

    sum(((d1[:, i, np.newaxis] - true_yield[i]) @
         (d1[:, i, np.newaxis] - true_yield[i]).T for i in range(d1.shape[1]))) / 14
    np.sqrt(_)
    # i = 0
    # (data[np.newaxis, i, :] - data[np.newaxis, i, :].mean()
    #  ).T @ (data[np.newaxis, i, :] - data[np.newaxis, i, :].mean())

    # @ (data[:, 1, np.newaxis] - data[:, 1, np.newaxis].mean()) / 14

    x1V = cov_matrix[0, 0]
    x2V = cov_matrix[1, 1]
    cov = cov_matrix[0, 1]

    BLUE_estimator(x1, x1V, x2, x2V, cov)
