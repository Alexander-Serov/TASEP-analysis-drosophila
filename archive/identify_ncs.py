

import copy
import logging
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import trange

from constants import (intensity_threshold, max_nc13_length_minutes,
                       max_reasonable_polymerase_number,
                       min_nc13_length_minutes, mins_per_frame, nc13_folder)
from Dataset import get_avg_on_regular_time_mesh
from save_series_plot import save_series_plot


def identify_ncs(data):
    """
    The new algorithm that works with Ben's data format.
    Calculate the average trace, then
    """
    datasets_len = data.dataset_id.max() + 1
    ncs_locations = pd.DataFrame(
        columns=['dataset_id', 'nc13_start', 'nc13_end', 'nc14_start'], index=range(0, datasets_len))

    # Output array
    data_out = copy.deepcopy(data)
    data_out['nc'] = np.nan

    # nc time limits dataframe
    iterables = pd.MultiIndex.from_product([range(datasets_len), range(11, 15)])
    nc_limits = pd.DataFrame(columns=['start', 'end'], index=iterables)

    # %% Identification of starts and ends of the ncs 13 and 14
    # TODO: Only 14 for the moment
    for dataset_id in trange(datasets_len, desc='Processing data sets'):
        # dataset_id = 10
        # start_nc13_frame = np.nan
        # end_nc13_frame = np.nan

        # Select dataset data
        dataset_data = data[(data.dataset_id == dataset_id)]

        # Calculate the average trace
        avg_trace = get_avg_on_regular_time_mesh(dataset_data, dt=0.7)
        # return
        # if dataset_id == 10:
        #     avg_trace.to_csv('temp.csv')

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

        # %% Now use the detected cycles to perform detection on the irregular mesh
        # Add as an extra column for the input data frame

        for nc in range(11, 15):
            nc_data = avg_trace[avg_trace.nc == nc]
            if not nc_data.empty:
                start = avg_trace[avg_trace.nc == nc].index[0]
                end = avg_trace[avg_trace.nc == nc].index[-1]
                nc_limits.loc[(dataset_id, nc), :] = [start, end]
        # print(nc_limits)

        # %% Assign ncs to original data

        for nc in range(11, 15):
            rows = ((data_out.time >= nc_limits.loc[(dataset_id, nc), 'start'])
                    & (data_out.time <= nc_limits.loc[(dataset_id, nc), 'end'])
                    & (data_out.dataset_id == dataset_id)
                    )
            data_out.loc[rows, 'nc'] = nc
    print(nc_limits)
    return data_out
    # print(list(data_out.nc.values))

    # # %% Plots
    # fig, ax = plt.subplots(1, 1, num=1, clear=True)
    # ax.plot(avg_trace.intensity)
    # xlims = plt.xlim()
    # ax.plot(xlims, [expr_threshold] * 2, 'r')
    #
    # fig, ax = plt.subplots(1, 1, num=2, clear=True)
    # ax.plot(avg_trace.nc)

    # %%

    # nc14_sum = ((avg_trace.nc >= 0) & (avg_trace.index >=0))
    # last_label = np.nan
    # for a in avg_trace:
    #     print(a)
    # print(avg_trace.index[2])
    # for i in trange(time_len):
    #     # print(thrsh_data.iloc[i].intensity)
    #     if avg_trace.loc[avg_trace.index[i], 'is_expressing'] == False:
    #         continue
    #     if avg_trace.loc[avg_trace.index[i - 1], 'is_expressing'] == False and not np.isnan(last_label):
    #         last_label += 1
    #     elif np.isnan(last_label) and float(avg_trace.index[i]) > 0:
    #         last_label = 14
    #     avg_trace.loc[avg_trace.index[i], 'nc'] = last_label
    # print(avg_trace)
    # print(list(avg_trace.nc.values))

    # start_nc14_frame = dataset_nc_data[dataset_nc_data.nc == 14].frame.iloc[0]
    # print(start_nc14_frame)


1
# %%
# # Manual nc locations corrections list for different data sets
# start_frame_corrections = {
#     '2014-10-13-hbBAC_NoPrim_B': 39,
#     '2014-10-14-hbBAC_NoPrim_B': 16,
#     '2014-10-30-hbBAC_NoPrim_A': 14,
#     '2014-07-08-kniBAC': 49,
#     '2014-07-09-kniBAC': 44,
#     '2014-03-01-SnaBAC_B': 41,
#     '2014-03-05-SnaBAC_NoPrim_B': 42,
#     '2014-04-06-SnaBAC_NoPrim_A': 41,
#     '2014-05-17-SnaBAC_NoShad_A': 62,
#     '2014-06-06-SnaBACA': 15,
#     '2014-06-10-SnaBAC_NoPrim_C': 30,
#     '2014-06-05-SnaBACA': 21,
# }
#
# start_nc13_frame_corrections = {
#     '2014-03-16-HbBAC_NoPrim_A': 8,
#     '2014-10-24-hbBAC_NoShad_C': np.nan,
#     '2014-08-20-kniBAC_NoShad_B': np.nan,
#     '2014-08-21-kniBAC_NoShad_A': np.nan,
#     '2014-05-27-SnaBAC_NoShad_A': np.nan,
#     '2014-05-28-SnaBAC_NoShad_A': np.nan,
#     '2014-06-05-SnaBACB': np.nan,
#     '2014-08-20-kniBAC_NoShad_A': np.nan,
#     '2014-08-22-kniBAC_NoShad_A': 17,
#     '2014-08-14-kniBAC_NoPrim_A': np.nan,
#     '2014-08-14-kniBAC_NoPrim_D': np.nan,
#     '2014-08-15-kniBAC_NoPrim_A': np.nan,
# }
# end_nc13_frame_corrections = {
#     '2014-03-16-HbBAC_NoPrim_A': 27,
#     '2014-10-24-hbBAC_NoShad_C': 4,
#     '2014-08-20-kniBAC_NoShad_B': 7,
#     '2014-08-21-kniBAC_NoShad_A': 4,
#     '2014-05-27-SnaBAC_NoShad_A': 6,
#     '2014-05-28-SnaBAC_NoShad_A': 5,
#     '2014-06-05-SnaBACB': 5,
#     '2014-08-20-kniBAC_NoShad_A': 35,
#     '2014-08-22-kniBAC_NoShad_A': 49,
#     '2014-08-14-kniBAC_NoPrim_A': np.nan,
#     '2014-08-14-kniBAC_NoPrim_D': np.nan,
#     '2014-08-15-kniBAC_NoPrim_A': 60,
# }
#
# # start_nc14_frame_corrections = {
# #     '2014-08-14-kniBAC_NoPrim_A': np.nan,
# # }
#
#
# # class NucCycLocationWarning(RuntimeWarning):
# #     NucCycLocationWarning
# #     pass
#
#
#
#
#
# def identify_ncs_madhavs_data(data):
#     """
#     Locate the start of nc13, nc14 and the end of nc13 to correctly measure the slopes
#     """
#
#     datasets_len = data.dataset_id.max() + 1
#     ncs_locations = pd.DataFrame(
#         columns=['dataset_id', 'nc13_start', 'nc13_end', 'nc14_start'], index=range(0, datasets_len))
#
#     skipped_polymerases_too_high = []
#
#     # %% Identification of starts and ends of the ncs 13 and 14
#     # dataset_id = 2
#     for dataset_id in trange(datasets_len):
#         # dataset_id = 29
#         [dataset_id]
#
#         start_nc13_frame = np.nan
#         end_nc13_frame = np.nan
#         start_nc14_frame = np.nan
#         bl_search_nc13 = True
#
#         # Select data
#         dataset_nc_data = data[(data.dataset_id == dataset_id)]
#
#         dataset_name = dataset_nc_data.dataset.iloc[0]
#         # print(dataset_id, dataset_name)
#
#         # Average
#         avg_dataset_nc_data = dataset_nc_data.groupby('Frame').mean()
#         # print(avg_dataset_nc_data)
#         # break
#
#         # Manual corrections for nc 13
#         if dataset_name in start_nc13_frame_corrections.keys():
#             start_nc13_frame = start_nc13_frame_corrections[dataset_name]
#             end_nc13_frame = end_nc13_frame_corrections[dataset_name]
#             bl_search_nc13 = False
#
#         # Manual corrections for first frame search
#         if dataset_name in start_frame_corrections.keys():
#             start_frame = start_frame_corrections[dataset_name]
#         else:
#             start_frame = avg_dataset_nc_data.index.min()
#
#         end_frame = avg_dataset_nc_data.index.max()
#
#         def has_data(frame):
#             return avg_dataset_nc_data[avg_dataset_nc_data.index == frame].count()[0] > 0
#
#         def check_nc_label(frame, var_name, expected_label):
#             if np.isfinite(frame):
#                 nc_label = dataset_nc_data.groupby('Frame').median().loc[frame].nc
#                 if nc_label != expected_label:
#                     logging.warning(
#                         f"Wrong nc label for {var_name} in dataset={dataset_id}: {nc_label}. Expected: {expected_label}")
#
#         # Search for the location of the ncs
#         tries = 5
#         if bl_search_nc13:
#             for i in range(tries):
#                 # Detect nc13 start frame as the start of first 3 sequential frames recorded
#                 for start_nc13_frame in range(start_frame, end_frame + 1):
#                     if has_data(start_nc13_frame) and has_data(start_nc13_frame + 1) and has_data(start_nc13_frame + 2):
#                         break
#                 start_nc13_frame
#
#                 check_nc_label(start_nc13_frame, 'nc13_start', 13)
#                 # # print(avg_dataset_nc_data)
#                 # nc_label = avg_dataset_nc_data.loc[start_nc13_frame].nc
#                 # if nc_label != 13:
#                 #     logging.warning(
#                 #         f"Wrong nc label for nc13_start in dataset={dataset_id}: {nc_label}")
#
#                 # Detect the end of the nc 13 by the gap in data
#                 for end_nc13_frame in range(start_nc13_frame, end_frame + 1):
#                     if not has_data(end_nc13_frame) and avg_dataset_nc_data[avg_dataset_nc_data.index < end_nc13_frame].count()[0] > 2:
#                         end_nc13_frame -= 1
#                         break
#                 end_nc13_frame
#
#                 check_nc_label(end_nc13_frame, 'nc13_end', 13)
#
#                 # If nc13 is too short, restart search at its end
#                 nc13_len_minutes = (end_nc13_frame - start_nc13_frame) * mins_per_frame
#                 nc13_len_minutes
#                 if nc13_len_minutes >= min_nc13_length_minutes:
#                     break
#                 else:
#                     start_frame = end_nc13_frame
#                     start_nc13_frame = np.nan
#         start_nc13_frame
#
#         # # If the next point is lower, choose it instead
#         avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index == start_nc13_frame + 1]
#         tries = 10
#         if not np.isnan(start_nc13_frame):
#             for i in range(tries):
#                 cur_polym = avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index ==
#                                                             start_nc13_frame].iloc[0]
#                 cur_polym
#
#                 next_frame = avg_dataset_nc_data.index[avg_dataset_nc_data.index > start_nc13_frame].min(
#                 )
#                 next_frame
#                 next_polym = avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index ==
#                                                              next_frame].iloc[0]
#                 next_polym
#
#                 if cur_polym >= next_polym:
#                     # Move to the next present frame
#                     start_nc13_frame = next_frame
#                 else:
#                     break
#         start_nc13_frame
#
#         # If the nc 13 is very long, it is nc 14 and not nc 13
#         nc13_len_minutes = (end_nc13_frame - start_nc13_frame) * mins_per_frame
#         if nc13_len_minutes > max_nc13_length_minutes:
#             start_nc14_frame = start_nc13_frame
#             end_nc13_frame = start_nc14_frame - 1
#             start_nc13_frame = np.nan
#         [start_nc13_frame, end_nc13_frame]
#
#         # Detect the start of the nc 14
#         if np.isfinite(end_nc13_frame):
#             for frame in range(end_nc13_frame + 1, end_frame + 1):
#                 if avg_dataset_nc_data[avg_dataset_nc_data.index == frame].count()[0] > 0 and avg_dataset_nc_data[avg_dataset_nc_data.index == frame + 1].count()[0] > 0 and avg_dataset_nc_data[avg_dataset_nc_data.index == frame + 2].count()[0] > 0:
#                     start_nc14_frame = frame
#                     break
#         else:
#             start_nc14_frame = np.nan
#
#         # # If in nc 14, the next point is lower, choose it
#         avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index == start_nc14_frame + 1]
#         tries = 10
#         if not np.isnan(start_nc14_frame):
#             for i in range(tries):
#                 cur_polym = avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index ==
#                                                             start_nc14_frame].iloc[0]
#                 cur_polym
#
#                 next_frame = avg_dataset_nc_data.index[avg_dataset_nc_data.index > start_nc14_frame].min(
#                 )
#                 next_frame
#                 next_polym = avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index ==
#                                                              next_frame].iloc[0]
#                 next_polym
#
#                 if cur_polym >= next_polym:
#                     # Move to the next present frame
#                     start_nc14_frame = next_frame
#                 else:
#                     break
#         start_nc14_frame
#
#         check_nc_label(start_nc14_frame, 'nc14_start', 14)
#
#         # If the dataset features polymerase numbers higher than 200, that must be a data error
#         # The data will be rechecked
#         if avg_dataset_nc_data.polymerases.max() >= max_reasonable_polymerase_number:
#             [start_nc13_frame, end_nc13_frame,
#                 start_nc14_frame] = np.nan * np.asarray([1, 1, 1])
#             skipped_polymerases_too_high.append(dataset_name)
#
#             # Save ncs locations
#         ncs_locations.loc[dataset_id] = [dataset_id,
#                                          start_nc13_frame, end_nc13_frame, start_nc14_frame]
#
#         # Plot
#         fit_interval = [0, 0]  # np.asarray(fit_frames) * mins_per_frame * 0
#         vbars = np.asarray([start_nc13_frame, end_nc13_frame, start_nc14_frame]) * mins_per_frame
#         filename = 'ncs_locations_%i.png' % (dataset_id)
#         filepath = os.path.join(nc13_folder, filename)
#         save_series_plot(x=avg_dataset_nc_data.time_min,
#                          y=avg_dataset_nc_data.polymerases, a=0, b=0, fit_interval=fit_interval, vbars=vbars, filename=filepath, dataset_name=dataset_name, dataset_id=dataset_id)
#         # return
#
#     # Report errors
#     if len(skipped_polymerases_too_high) > 0:
#         # warnings.warn
#         logging.warning("%i datasets skipped due to polymerase values >= %.0f:\n" %
#                         (len(skipped_polymerases_too_high), max_reasonable_polymerase_number))
#         print(skipped_polymerases_too_high)
#
#     return ncs_locations
