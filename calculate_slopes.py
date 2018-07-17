

import numpy as np
import pandas as pd

from constants import slope_length_frames


def calculate_slopes(data, ncs_locations, trace_id):
    """
    Calculate nc13 and nc14 slope for the given trace_id
    """
    dataset_id = data.dataset_id.iloc[0]
    gene_id = data.gene_id.iloc[0]
    construct_id = data.construct_id.iloc[0]
    dataset_id

    def get_slope(start_frame):
        end_slope_frame = start_frame + slope_length_frames - 1
        cur_data = data[(data.Frame >= start_frame) & (data.Frame <= end_slope_frame)]
        if cur_data.count().Frame <= 1 or np.isnan(start_frame):
            slope = np.nan
        else:
            coefs = np.polyfit(cur_data.time_min, cur_data.polymerases, 1)
            slope = coefs[0]
        return slope

    def get_mean_max(start_frame, end_frame):
        cur_data = data[(data.Frame >= start_frame) & (data.Frame <= end_frame)]
        if cur_data.count().Frame < 3 or np.isnan(start_frame) or np.isnan(end_frame):
            max_polymerase = np.nan
        else:
            max_frame = cur_data[cur_data.polymerases ==
                                 cur_data.polymerases.max()].Frame.iloc[0]
            # Get the average in 3 frames around the maximum. If there are no frames around, skip
            data_around_max = cur_data[(cur_data.Frame >= max_frame - 1) &
                                       (cur_data.Frame <= max_frame + 1)]
            if data_around_max.count().Frame < 3:
                max_polymerase = np.nan
            else:
                max_polymerase = data_around_max.mean().polymerases
        return max_polymerase

    # nc 13 slope
    start_nc13_frame = ncs_locations.iloc[dataset_id].nc13_start
    slope_nc13 = get_slope(start_nc13_frame)

    # nc 13 max value
    end_nc13_frame = ncs_locations.iloc[dataset_id].nc13_end
    max_polymerase_nc13 = get_mean_max(start_nc13_frame, end_nc13_frame)

    # nc 14 slope
    start_nc14_frame = ncs_locations.iloc[dataset_id].nc14_start
    slope_nc14 = get_slope(start_nc14_frame)

    # nc 14 max value
    end_nc14_frame = data[(data.trace_id == trace_id)].max().Frame
    max_polymerase_nc14 = get_mean_max(start_nc14_frame, end_nc14_frame)

    slope_series = pd.Series({'dataset_id': dataset_id, 'gene_id':
                              gene_id, 'construct_id': construct_id, 'slope_nc13': slope_nc13, 'slope_nc14': slope_nc14, 'max_nc13': max_polymerase_nc13, 'max_nc14': max_polymerase_nc14})
    return slope_series
