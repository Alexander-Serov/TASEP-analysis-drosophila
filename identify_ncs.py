

import copy
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import trange

from constants import (default_intensity_threshold, dt_new,
                       intensity_thresholds, min_nc13_length_minutes,
                       mins_per_frame, nc13_folder)
from Dataset import get_avg_on_regular_time_mesh

# In some data sets, nc14 was absent and only nc13 was present. Needed to indentify ncs differently
# By default the last observed nc in nc14. Add elements to the dictionary if manual adjustements are necessary
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
# 2: 15,
#             6: 16,
#             8: 15,
#             10: 15,
# 41: 15, 45: 15, 56: 15, 59: 15, 66: 15, 67: 15}

# dataset_ids_with_no_nc14 = (20, 21, 22)


def identify_ncs(data):
    """
    The new algorithm that works with Ben's data format.
    NC identification is implemented as first identifying continuous expression regions over the given threshold. Then they are continuously numbered with a cumsum.
    Assume then that the last observed nc is nc14 and shift everything accordingly.
    """
    datasets_len = data.dataset_id.max() + 1

    # Prepare an output df
    data_out = copy.deepcopy(data)
    data_out['nc'] = np.nan

    # Create nc time limits dataframe
    iterables = pd.MultiIndex.from_product([range(datasets_len), range(11, 15)])
    nc_limits = pd.DataFrame(columns=['Tstart', 'Tend'], index=iterables)

    for dataset_id in trange(datasets_len, desc='Processing data sets'):
        # dataset_id = 13
        dataset_data = data[(data.dataset_id == dataset_id)]
        # print(dataset_id, dataset_data)
        avg_trace, _ = get_avg_on_regular_time_mesh(dataset_data, dt=dt_new)

        # %% Assign nc numbers assuming that the first one after time t > 0 is nc14
        intensity_threshold = intensity_thresholds.get(dataset_id, default_intensity_threshold)
        avg_trace['is_expressing'] = avg_trace >= intensity_threshold
        # Find nc Tstarts
        avg_trace['nc'] = 1 * (avg_trace.is_expressing &
                               np.logical_not(avg_trace.is_expressing.shift(1)))
        # Take care of the first time point
        avg_trace.loc[avg_trace.index[0], 'nc'] = 1 * \
            avg_trace.loc[avg_trace.index[0], 'is_expressing']

        # Convert group Tstarts to nc numbers
        avg_trace['nc'] = avg_trace['nc'].cumsum() * avg_trace.is_expressing
        # print(intensity_threshold)
        # print(avg_trace)

        # Assume nc14 is the last nc present
        # try:
        # avg_trace[avg_trace.is_expressing & (avg_trace.index >= 0)].nc.iloc[0]
        if dataset_id in last_ncs:
            last_nc = last_ncs[dataset_id]
        else:
            last_nc = 14
        nc14_offset = avg_trace.nc.max()
        avg_trace['nc'] += (last_nc - nc14_offset)
        # print(avg_trace)

        # plt.figure(num=2, clear=True)
        # # plt.scatter(dataset_data.time, dataset_data.intensity)
        # plt.plot(avg_trace.index, avg_trace.intensity, '-o')
        # plt.plot([0, avg_trace.index.max()], [intensity_threshold] * 2, 'g')
        # plt.show()

        # break

        # Drop all ncs that are shorter than 3 time points
        # avg_trace

        # print(avg_trace)
        # break
        # except:
        # if dataset_id not in dataset_ids_with_no_nc14:
        #     logging.warning(
        #         f'nc14 not found in dataset {dataset_id}. Will use nc13. Please check')
        # nc13_offset = avg_trace[avg_trace.is_expressing &
        #                         (avg_trace.index <= 0)].nc.tail(1).values
        # # print(nc13_offset)
        # avg_trace['nc'] += (13 - nc13_offset)
        avg_trace.loc[~avg_trace.is_expressing, 'nc'] = np.nan
        # print(avg_trace)

        # %% Collect the nc time intervals in a separate table
        # Use the detected nc to label the data on the irregular mesh
        # Index key is a tuple
        # Add as an extra column for the input data frame
        for nc in range(11, 15):
            nc_data = avg_trace[avg_trace.nc == nc]
            if not nc_data.empty:
                Tstart = avg_trace[avg_trace.nc == nc].index[0]
                Tend = avg_trace[avg_trace.nc == nc].index[-1]
                nc_limits.loc[(dataset_id, nc), ['Tstart', 'Tend']] = [Tstart, Tend]

        # %% Use the table to label original non-averaged data
        for nc in range(11, 15):
            rows = ((data_out.time >= nc_limits.loc[(dataset_id, nc), 'Tstart'])
                    & (data_out.time <= nc_limits.loc[(dataset_id, nc), 'Tend'])
                    & (data_out.dataset_id == dataset_id)
                    )
            data_out.loc[rows, 'nc'] = nc
    # print(nc_limits)
    return data_out, nc_limits
