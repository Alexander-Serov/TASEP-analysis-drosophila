

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from constants import mins_per_frame, output_slopes_folder


def find_best_fit(data):
    """
    Identify the location of the initial polymerase injection slope by finding the lowest residue fit of a given length to the given data points. The same number of frames is always taken into account, but not necessarily the same time interval (if some frames are missing)
    """

    avg_slope_length_frames = 8
    # %% Initialize
    start_frame = data.Frame.min()
    end_frame = data.Frame.max()
    norm_residual_best = np.inf

    # Check if contains more frames than expected for fitting
    if end_frame - start_frame + 1 >= avg_slope_length_frames:
        # Parse all possible fit start frames_len

        # print(start_frame, end_frame)

        for fit_start_frame in range(start_frame, end_frame - avg_slope_length_frames + 2):
            fit_end_frame = fit_start_frame + avg_slope_length_frames - 1
            data_slice = data[(data.Frame >= fit_start_frame) & (data.Frame <= fit_end_frame)]
            points_count = data_slice.count()[0]

            # The data slice must contain at least two different frames
            if data_slice.Frame.nunique() <= 2:
                continue

            # print(data_slice)
            coefs, residual = np.polyfit(
                data_slice.time_min, data_slice.polymerases, 1, full=True)[0:2]
            norm_residual = residual[0] / points_count
            # print([fit_start_frame, fit_end_frame])
            # print(coefs, residual)
            # print([coefs[0], -coefs[1] / coefs[0], norm_residual])

            # Recalculate coefs
            intersect = -coefs[1] / coefs[0]
            slope = coefs[0]

            # Compare residual to the best residual so far given positive slope
            if norm_residual < norm_residual_best and slope > 0 and intersect > 0:
                norm_residual_best = norm_residual
                intersect_best = intersect
                slope_best = slope
                fit_frames = [fit_start_frame, fit_end_frame]

        # Calculate the approximate length of the initial slope
        approx_slope_length_min = end_frame * mins_per_frame - intersect_best

    elif end_frame - start_frame + 1 > 2:
        # If few points, but > 2 available, perform the only fit possible
        points_count = data.count()[0]
        coefs, residual = np.polyfit(
            data.time_min, data.polymerases, 1, full=True)[0:2]
        norm_residual = residual[0] / points_count

        # Compare residual to the best residual so far
        if norm_residual < norm_residual_best:
            norm_residual_best = norm_residual
            intersect_best = -coefs[1] / coefs[0]
            slope_best = coefs[0]
            fit_frames = [start_frame, end_frame]

    # Calculate the approximate length of the initial slope
    approx_slope_length_min = end_frame * mins_per_frame - intersect_best

    return [intersect_best, slope_best, approx_slope_length_min, fit_frames, norm_residual_best]


# def save_series_plot(x, y, a, b, hold=False):
#
#     fig = plt.figure(1)
#     if not hold:
#         fig.clf()
#     ax = plt.gca()
#
#     # Plot trace
#     ax.plot(x, y)
#
#     # Plot fit
#     print([a, b])
#     x_fit = np.asarray([np.min(x), x[y == np.max(y)]])
#     y_fit = x_fit * a + b
#     ax.plot(x_fit, y_fit, 'r')
#
#     plt.show
