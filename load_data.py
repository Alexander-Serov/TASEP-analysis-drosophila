

try:
    has_run
except NameError:
    %matplotlib tk
    %load_ext autoreload
    %autoreload 2
    has_run = 1
else:
    print("Graphic interface NOT re-initialized")


from find_best_fit import find_best_fit
from identify_shift_and_slope import identify_shift_and_slope
import numpy as np
import os
import pandas as pd
from reinit_folder import reinit_folder
from save_series_plot import save_series_plot
from tqdm import trange
from constants import matlab_csv_filepath, expected_slope_length_frames, mins_per_frame, output_slopes_folder, detect_nc13_leftovers_interval_frames, nc13_folder, min_nc13_length_minutes, max_nc13_length_minutes


# %% Initialize
reinit_folder(output_slopes_folder)
reinit_folder(nc13_folder)


# %% Import
data = pd.read_csv(matlab_csv_filepath, sep=';')
print(data.columns.values)

# %% Analyze
traces_len = data.trace_id.max() + 1
datasets_len = data.dataset_id.max() + 1
intersects = np.ones([traces_len, 1]) * np.nan
slopes = np.ones([traces_len, 1]) * np.nan

# # Calculate number of traces overall in each nc
# nc13_total_count, nc14_total_count = [data[data.nc == nc].trace_id.nunique() for nc in[13, 14]]
# print([nc13_total_count, nc14_total_count])


# %% Identification of starts and ends of the ncs 13 and 14
# nc = 13
dataset_id = 2
for dataset_id in trange(datasets_len):
    print([dataset_id, nc])

    start_nc13_frame = np.nan
    end_nc13_frame = np.nan
    start_nc14_frame = np.nan
    bl_search_nc13 = True

    # Select data
    # dataset = data[data.dataset_id == dataset_id]
    dataset_nc_data = data[(data.dataset_id == dataset_id)]

    # # If no nc 13, skip
    # if dataset_nc_data.count()[0] == 0:
    #     continue

    dataset_name = dataset_nc_data.dataset.iloc[0]
    dataset_name

    # Average
    avg_dataset_nc_data = dataset_nc_data.groupby('Frame').mean()
    avg_dataset_nc_data
    start_frame = avg_dataset_nc_data.index.min()
    end_frame = avg_dataset_nc_data.index.max()

    # Manual corrections
    if dataset_id == 2:
        start_nc13_frame = 8
        end_nc13_frame = 27
        bl_search_nc13 = False
    if dataset_id == 31:
        start_frame = 39
    elif dataset_id == 32:
        start_frame = 16
    elif dataset_id == 36:
        start_nc13_frame = np.nan
        end_nc13_frame = 4
        bl_search_nc13 = False
    elif dataset_id == 38:
        start_frame = 14
    if dataset_id == 44:
        start_frame = 49
    elif dataset_id == 45:
        start_frame = 44
    elif dataset_id == 52:
        start_nc13_frame = np.nan
        end_nc13_frame = 7
        bl_search_nc13 = False
    elif dataset_id == 53:
        start_nc13_frame = np.nan
        end_nc13_frame = 4
        bl_search_nc13 = False
    elif dataset_id == 56:
        start_frame = 41
    elif dataset_id == 59:
        start_frame = 42
    elif dataset_id == 62:
        start_frame = 41
    elif dataset_id == 64:
        start_frame = 62
    elif dataset_id == 65:
        start_nc13_frame = np.nan
        end_nc13_frame = 6
        bl_search_nc13 = False
    elif dataset_id == 66:
        start_nc13_frame = np.nan
        end_nc13_frame = 5
        bl_search_nc13 = False
    elif dataset_id == 67:
        start_frame += 20
    elif dataset_id == 68:
        start_nc13_frame = np.nan
        end_nc13_frame = 5
        bl_search_nc13 = False
    elif dataset_id == 69:
        start_frame = 15
    elif dataset_id == 76:
        start_frame = 30
    [start_frame, end_frame]

    tries = 5
    if bl_search_nc13:
        for i in range(tries):
            # Detect nc13 start frame
            for start_nc13_frame in range(start_frame, end_frame + 1):
                if avg_dataset_nc_data[avg_dataset_nc_data.index == start_nc13_frame].count()[0] > 0 and avg_dataset_nc_data[avg_dataset_nc_data.index == start_nc13_frame + 1].count()[0] > 0 and avg_dataset_nc_data[avg_dataset_nc_data.index == start_nc13_frame + 2].count()[0] > 0:
                    break
            start_nc13_frame

            # Detect the end of the nc 13 by the gap in data
            for end_nc13_frame in range(start_nc13_frame, end_frame + 1):
                if avg_dataset_nc_data[avg_dataset_nc_data.index == end_nc13_frame].count()[0] == 0 and avg_dataset_nc_data[avg_dataset_nc_data.index < end_nc13_frame].count()[0] > 2:
                    end_nc13_frame -= 1
                    break
            end_nc13_frame

            # If nc13 is too short, restart search at its end
            nc13_len_minutes = (end_nc13_frame - start_nc13_frame) * mins_per_frame
            nc13_len_minutes
            if nc13_len_minutes >= min_nc13_length_minutes:
                break
            else:
                start_frame = end_nc13_frame
                start_nc13_frame = np.nan
    start_nc13_frame

    # # If the next point is lower, choose it instead
    avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index == start_nc13_frame + 1]
    tries = 10
    if not np.isnan(start_nc13_frame):
        for i in range(tries):
            cur_polym = avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index ==
                                                        start_nc13_frame].iloc[0]
            cur_polym

            next_frame = avg_dataset_nc_data.index[avg_dataset_nc_data.index > start_nc13_frame].min(
            )
            next_frame
            next_polym = avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index ==
                                                         next_frame].iloc[0]
            next_polym

            if cur_polym >= next_polym:
                # Move to the next present frame
                start_nc13_frame = next_frame
            else:
                break
    start_nc13_frame

    # If the nc 13 is very long, it is nc 14 and not nc 13
    nc13_len_minutes = (end_nc13_frame - start_nc13_frame) * mins_per_frame
    if nc13_len_minutes > max_nc13_length_minutes:
        start_nc14_frame = start_nc13_frame
        end_nc13_frame = start_nc14_frame - 1
        start_nc13_frame = np.nan
    [start_nc13_frame, end_nc13_frame]

    # Detect the start of the nc 14
    for frame in range(end_nc13_frame + 1, end_frame + 1):
        if avg_dataset_nc_data[avg_dataset_nc_data.index == frame].count()[0] > 0 and avg_dataset_nc_data[avg_dataset_nc_data.index == frame + 1].count()[0] > 0 and avg_dataset_nc_data[avg_dataset_nc_data.index == frame + 2].count()[0] > 0:
            start_nc14_frame = frame
            break

    # # If in nc 14, the next point is lower, choose it
    avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index == start_nc14_frame + 1]
    tries = 10
    if not np.isnan(start_nc14_frame):
        for i in range(tries):
            cur_polym = avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index ==
                                                        start_nc14_frame].iloc[0]
            cur_polym

            next_frame = avg_dataset_nc_data.index[avg_dataset_nc_data.index > start_nc14_frame].min(
            )
            next_frame
            next_polym = avg_dataset_nc_data.polymerases[avg_dataset_nc_data.index ==
                                                         next_frame].iloc[0]
            next_polym

            if cur_polym >= next_polym:
                # Move to the next present frame
                start_nc14_frame = next_frame
            else:
                break
    start_nc14_frame

    # # Manual corrections
    # if dataset_id == 26:
    #     start_nc14_frame += 1
    # elif dataset_id == 29:
    #     start_nc14_frame += 1
    # elif dataset_id == 54:
    #     start_nc14_frame += 4
    # elif dataset_id == 69:
    #     start_nc14_frame += 1

    # elif dataset_id == 55:
    #     start_nc14_frame += 4
    # elif dataset_id == 74:
    #     start_nc14_frame += 2
    start_nc14_frame

    # Plot
    fit_interval = [0, 0]  # np.asarray(fit_frames) * mins_per_frame * 0
    vbars = np.asarray([start_nc13_frame, end_nc13_frame, start_nc14_frame]) * mins_per_frame
    filename = 'slopes_dataset_%i_nc_%i.png' % (dataset_id, nc)
    filepath = os.path.join(nc13_folder, filename)
    save_series_plot(x=avg_dataset_nc_data.time_min,
                     y=avg_dataset_nc_data.polymerases, a=0, b=0, fit_interval=fit_interval, vbars=vbars, filename=filepath, dataset_name=dataset_name)


# %% nc 14 analysis
analyzed_count = 0
processed_successfully = 0
skipped_short_slope = 0
skipped_too_short = 0
skipped_negative_slope = 0
# skipped_starts_with_the_maximum = 0
skipped_discontinuous = 0

short_trace_lengths = []
short_slope_lengths = []
negative_slopes = []

for dataset_id in trange(datasets_len):
    # dataset_id = 72
    nc = 14
    print([dataset_id, nc])
    # Select data
    dataset = data[data.dataset_id == dataset_id]
    dataset_nc_data = dataset[dataset.nc == nc]
    dataset.columns.values
    dataset_name = dataset.dataset.iloc[0]

    # Skip if only nc13 is present
    dataset_nc_data.count()[0]
    if dataset_nc_data.count()[0] < 1:
        continue

    dataset_nc_data

    # Average the trace dataset in nc14
    avg_dataset_nc_data = dataset_nc_data.groupby('Frame').mean()

    # Skip if less than 2 frames are present
    if len(avg_dataset_nc_data) < 2:
        continue
    avg_dataset_nc_data

    # Detect if there are leftovers from nc13 by looking at gaps in frames
    start_frame = avg_dataset_nc_data.index.min()
    start_frame
    expected_end_nc13_frame = start_frame + detect_nc13_leftovers_interval_frames - 1
    for end_nc13_frame in range(expected_end_nc13_frame, start_frame - 1, -1):
        if avg_dataset_nc_data[avg_dataset_nc_data.index == end_nc13_frame].count()[0] < 1:
            break
    end_nc13_frame

    # Cut what's left of the nc 13
    avg_dataset_nc_data = avg_dataset_nc_data[avg_dataset_nc_data.index > end_nc13_frame]
    avg_dataset_nc_data

    max_frame = avg_dataset_nc_data[avg_dataset_nc_data.polymerases ==
                                    avg_dataset_nc_data.polymerases.max()].index[0]
    max_frame

    # Manual correction for some data sets
    if dataset_id == 21:
        max_frame = start_frame + 12
    elif dataset_id == 72:
        max_frame = start_frame + 12

        # # If the maximum is at the first frame, the nc13 has not been properly cut off.
        # # Repeat search without the first point
        # # Drop by index
        #     if max_frame == avg_dataset_nc_data.index[0]:
        #         avg_dataset_nc_data = avg_dataset_nc_data.drop(max_frame)
        #     else:
        #         bl_max_found = True

        # earlier_frames = avg_dataset_nc_data.iloc[avg_dataset_nc_data.index < max_frame]
        # earlier_frames
        # min_frame = earlier_frames[earlier_frames.polymerases ==
        #                            earlier_frames.polymerases.min()].index[0]
        # earlier_frames[earlier_frames.polymerases == earlier_frames.polymerases.min()]

        # Start search at the end of nc13
    min_frame = end_nc13_frame

    # print(min_frame, max_frame)

    # % Fit the average trace between these frames with a sliding window of fixed length and find the best fit
    potential_slope_data = dataset_nc_data[(dataset_nc_data.Frame >= min_frame) & (
        dataset_nc_data.Frame <= max_frame)][['Frame', 'polymerases', 'dataset_id', 'nc', 'time_min']]

    intersect, slope, approx_slope_length_min, fit_frames, _ = find_best_fit(
        potential_slope_data)
    print(approx_slope_length_min)
    # print([intersect, slope, approx_slope_length_min])
    # potential_slope_data

    # print(avg_dataset_nc_data.index)

    # Plot with fit
    fit_interval = np.asarray(fit_frames) * mins_per_frame
    filename = 'slopes_dataset_%i_nc_%i.png' % (dataset_id, nc)
    save_series_plot(x=avg_dataset_nc_data.time_min,
                     y=avg_dataset_nc_data.polymerases, a=slope, b=-intersect * slope, fit_interval=fit_interval, filename=filename, dataset_name=dataset_name)

# %%
1
# %%
1
# # for trace_id in trange(1 + 0 * traces_len):
# analyzed_count += 1
# cur_trace = data[data.trace_id == trace_id]
#
# frames_len = cur_trace.count().Frame
#
# # print(frames_len)
# # A trace must contain at least `expected_slope_length_frames` frames
# if frames_len < expected_slope_length_frames:
#     skipped_too_short += 1
#     short_trace_lengths.append(frames_len)
#     continue
#
# start_frame = cur_trace.Frame.iloc[0]
# # Identify the first maximum fluo for the given trace
# max_frame = cur_trace.Frame[cur_trace.polymerases == cur_trace.polymerases.max()].iloc[0]
#
# # There must be at least `expected_slope_length_frames` before the maximum
# earlier_frames = cur_trace[cur_trace.Frame <= max_frame]
# if earlier_frames.count().Frame < expected_slope_length_frames:
#     skipped_short_slope += 1
#     short_slope_lengths.append(earlier_frames.count().Frame)
#     continue
#
# # Get the possible locations for the start of the fixed-frame count fit
# earlier_frames = earlier_frames.iloc[0:(-expected_slope_length_frames + 1)]
#
# # Find the last minimum among these locations
# min_frame = earlier_frames[earlier_frames.polymerases ==
#                            earlier_frames.polymerases.min()].Frame.iloc[-1]
#
# # Extract frames with the potential slope between the last minimum and the maximum
# potential_slope_frames = cur_trace[(cur_trace.Frame >= min_frame)
#                                    & (cur_trace.Frame <= max_frame)]
#
# slope, intersect, residual = identify_shift_and_slope(potential_slope_frames.time_min.values,
#                                                       potential_slope_frames.polymerases.values, potential_slope_frames.nc.iloc[0])
# # print([slope, intersect, residual])
#
# # Drop if slope is negative
# if slope < 0:
#     skipped_negative_slope += 1
#     negative_slopes.append(slope)
#     continue
#
# # Store
# slopes[trace_id] = slope
# intersects[trace_id] = intersect
# processed_successfully += 1

# # Slope detection statistics
# print('%i traces analyzed' % analyzed_count)
# print('%i traces successfully processed' % processed_successfully)
# print('%i traces skipped because too short. Median length: %.i' %
#       (skipped_too_short, np.median(short_trace_lengths)))
# # print('%i traces skipped because started with the maximum' % skipped_starts_with_the_maximum)
# print('%i traces skipped because of short slope. Median slope length: %.i' %
#       (skipped_short_slope, np.median(short_slope_lengths)))
# print('%i traces skipped because negative slope. Median slope: %.2f' %
#       (skipped_negative_slope, np.median(negative_slopes)))
# print('%i traces skipped because discontinuous' % skipped_discontinuous)


# %%
print(slopes.tolist())
print(intersects.tolist())
# print(data[data.trace_id == 56])
print(np.nanmedian(data.APpos))
