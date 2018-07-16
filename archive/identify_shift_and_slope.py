import numpy as np
from constants import expected_slope_length_frames


def identify_shift_and_slope(time_mins, pols, nc):
    """
    Identify the location of the initial polymerase injection slope by finding the lowest residue fit of a given length to the given data points. The same number of frames is always taken into account, but not necessarily the same time interval (if some frames are missing)
    """

    # %% Initialize
    # slope_best = []
    # intersect_best = []
    residual_best = np.inf
    frames_len = len(time_mins)

    # Go through the possible slope locations
    # The number of points is fixed, but not necessarily the time interval
    for start_frame in range(frames_len - expected_slope_length_frames + 1):
        # Find the best fit (least squared difference) of a fixed length to the data
        coefs, residual = np.polyfit(time_mins[start_frame:(start_frame + expected_slope_length_frames)],
                                     pols[start_frame:(
                                         start_frame + expected_slope_length_frames)],
                                     1, full=True)[0:2]

        # Compare residual to the best residual so far
        if residual[0] < residual_best:
            residual_best = residual[0]
            slope_best = coefs[0]
            intersect_best = -coefs[1] / coefs[0]

    return [slope_best, intersect_best, residual_best]
