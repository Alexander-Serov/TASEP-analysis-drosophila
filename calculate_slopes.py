

import numpy as np
import pandas as pd

from constants import slope_length_frames


def calculate_slopes(data, ncs_locations, trace_id):
    """
    Calculate nc13 and nc14 slope for the given trace_id
    """
    dataset_id = data[data.trace_id == trace_id].dataset_id.iloc[0]
    gene_id = data[data.trace_id == trace_id].gene_id.iloc[0]
    construct_id = data[data.trace_id == trace_id].construct_id.iloc[0]
    dataset_id

    # nc 13
    start_nc13_frame = ncs_locations.iloc[dataset_id].nc13_start
    start_nc13_frame
    end_nc13_slope_frame = start_nc13_frame + slope_length_frames - 1
    end_nc13_slope_frame
    cur_nc13_data = data[(data.trace_id == trace_id) & (
        data.Frame >= start_nc13_frame) & (data.Frame <= end_nc13_slope_frame)]
    cur_nc13_data
    not cur_nc13_data.count().Frame or np.isnan(start_nc13_frame)
    if cur_nc13_data.count().Frame <= 1 or np.isnan(start_nc13_frame):
        slope_nc13 = np.nan
    else:
        coefs = np.polyfit(cur_nc13_data.time_min, cur_nc13_data.polymerases, 1)
        slope_nc13 = coefs[0]

    # nc 14
    start_nc14_frame = ncs_locations.iloc[dataset_id].nc14_start
    start_nc14_frame
    end_nc14_slope_frame = start_nc14_frame + slope_length_frames - 1
    end_nc14_slope_frame
    cur_nc14_data = data[(data.trace_id == trace_id) & (
        data.Frame >= start_nc14_frame) & (data.Frame <= end_nc14_slope_frame)]
    cur_nc14_data

    if cur_nc14_data.count().Frame <= 1 or np.isnan(start_nc14_frame):
        slope_nc14 = np.nan
    else:
        coefs = np.polyfit(cur_nc14_data.time_min, cur_nc14_data.polymerases, 1)
        slope_nc14 = coefs[0]
    [slope_nc13, slope_nc14]

    slope_series = pd.Series({'dataset_id': dataset_id, 'gene_id':
                              gene_id, 'construct_id': construct_id, 'slope_nc13': slope_nc13, 'slope_nc14': slope_nc14})
    return slope_series
