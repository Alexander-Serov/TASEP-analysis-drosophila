

import copy

import numpy as np
import pandas as pd
from tqdm import trange

from constants import default_intensity_threshold, dt_new, intensity_thresholds
from Dataset import get_avg_on_regular_time_mesh


def identify_ncs(data):
    """
    Identify ncs in data traces by first thresholding and then numbering the ncs above the threshold.
    See README form more details.
    """
    # The sequential number of the last observed nc.
    # By default the last expected nc is nc14.
    # Add elements to the following dictionary if manual adjustements are necessary
    last_ncs = {
        1: 15,
        6: 14,
        54: 13,
        55: 13,
        56: 13,
        57: 13,
        58: 14,
        59: 14,
        60: 12,
        67: 13,
        68: 13,
        70: 13,
    }

    datasets_len = data.dataset_id.max() + 1
    ncs = range(11, 15)

    # Prepare output data frames
    data_out = copy.deepcopy(data)
    data_out['nc'] = np.nan
    iterables = pd.MultiIndex.from_product(
        [range(datasets_len), ncs], names=['dataset_id', 'nc'])
    nc_limits = pd.DataFrame(columns=['Tstart', 'Tend'], index=iterables)

    for dataset_id in trange(datasets_len, desc='Processing data sets'):
        # dataset_id = 54
        dataset_data = data[(data.dataset_id == dataset_id)]

        # Calculate average trace
        avg_trace, _ = get_avg_on_regular_time_mesh(dataset_data, dt=dt_new)

        # Threshold the average trace
        intensity_threshold = intensity_thresholds.get(dataset_id, default_intensity_threshold)
        avg_trace['is_expressing'] = avg_trace >= intensity_threshold

        # Detect the start of an nc as a sequence of 2 points, where the first one is below the threshold and the second one is above
        # Assign 1 to all such nc starts
        avg_trace['nc'] = 1 * (avg_trace.is_expressing &
                               np.logical_not(avg_trace.is_expressing.shift(1)))
        # Assign 1 also to the first time point
        avg_trace.loc[avg_trace.index[0], 'nc'] = 1 * \
            avg_trace.loc[avg_trace.index[0], 'is_expressing']

        # Number the cycles by performing a cumsum
        avg_trace['nc'] = avg_trace['nc'].cumsum() * avg_trace.is_expressing

        # Shift the numbers based on the last nc number
        if dataset_id in last_ncs:
            last_nc = last_ncs[dataset_id]
        else:
            last_nc = 14
        nc14_offset = avg_trace.nc.max()
        avg_trace['nc'] += (last_nc - nc14_offset)

        # If a frame is not expressing, set to nan
        avg_trace.loc[~avg_trace.is_expressing, 'nc'] = np.nan

        # %% Collect the nc time intervals in a separate table
        for nc in ncs:
            nc_data = avg_trace[avg_trace.nc == nc]
            if not nc_data.empty:
                Tstart = avg_trace[avg_trace.nc == nc].index[0]
                Tend = avg_trace[avg_trace.nc == nc].index[-1]
                nc_limits.loc[(dataset_id, nc), ['Tstart', 'Tend']] = [Tstart, Tend]

        # %% Use the table to label original non-averaged data set
        for nc in ncs:
            indices = ((data_out.time >= nc_limits.loc[(dataset_id, nc), 'Tstart'])
                       & (data_out.time <= nc_limits.loc[(dataset_id, nc), 'Tend'])
                       & (data_out.dataset_id == dataset_id)
                       )
            # Add as an extra column for the input data frame
            data_out.loc[indices, 'nc'] = nc
    return data_out, nc_limits
