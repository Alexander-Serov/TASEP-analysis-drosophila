

import copy

import numpy as np
import pandas as pd
from tqdm import trange

from constants import (dt_new, intensity_threshold, min_nc13_length_minutes,
                       mins_per_frame, nc13_folder)
from Dataset import get_avg_on_regular_time_mesh


def identify_ncs(data):
    """
    The new algorithm that works with Ben's data format.
    NC identification is implemented as first identifying continuous expression regions over the given threshold. Then they are continuously numbered with a cumsum. Next the first one for t>0 is considered to be nc 14, and all numbers are shifted accrodingly
    """
    datasets_len = data.dataset_id.max() + 1

    # Prepare an output df
    data_out = copy.deepcopy(data)
    data_out['nc'] = np.nan

    # Create nc time limits dataframe
    iterables = pd.MultiIndex.from_product([range(datasets_len), range(11, 15)])
    nc_limits = pd.DataFrame(columns=['start', 'end'], index=iterables)

    for dataset_id in trange(datasets_len, desc='Processing data sets'):
        dataset_data = data[(data.dataset_id == dataset_id)]
        avg_trace = get_avg_on_regular_time_mesh(dataset_data, dt=dt_new)

        # %% Assign nc numbers assuming that the first one after time t > 0 is nc14
        avg_trace['is_expressing'] = avg_trace >= intensity_threshold
        # Find nc starts
        avg_trace['nc'] = 1 * (avg_trace.is_expressing &
                               np.logical_not(avg_trace.is_expressing.shift(1)))
        # Take care of the first time point
        avg_trace.loc[avg_trace.index[0], 'nc'] = 1 * \
            avg_trace.loc[avg_trace.index[0], 'is_expressing']

        # Convert group starts to nc numbers
        avg_trace['nc'] = avg_trace['nc'].cumsum() * avg_trace.is_expressing
        nc14_offset = avg_trace[avg_trace.is_expressing & (avg_trace.index >= 0)].nc.iloc[0]
        avg_trace['nc'] += (14 - nc14_offset)
        avg_trace.loc[~avg_trace.is_expressing, 'nc'] = np.nan

        # %% Collect the nc time intervals in a separate table
        Use the detected nc to label the data on the irregular mesh
        # Index key is a tuple
        # Add as an extra column for the input data frame
        for nc in range(11, 15):
            nc_data = avg_trace[avg_trace.nc == nc]
            if not nc_data.empty:
                start = avg_trace[avg_trace.nc == nc].index[0]
                end = avg_trace[avg_trace.nc == nc].index[-1]
                nc_limits.loc[(dataset_id, nc), ['start', 'end']] = [start, end]

        # %% Use the table to label original non-averaged data
        for nc in range(11, 15):
            rows = ((data_out.time >= nc_limits.loc[(dataset_id, nc), 'start'])
                    & (data_out.time <= nc_limits.loc[(dataset_id, nc), 'end'])
                    & (data_out.dataset_id == dataset_id)
                    )
            data_out.loc[rows, 'nc'] = nc
    return data_out
