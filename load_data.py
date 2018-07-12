

from identify_shift_and_slope import identify_shift_and_slope
import numpy as np
import pandas as pd
from tqdm import trange
from constants import matlab_csv_filepath, expected_slope_length_frames

try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = 1
else:
    print("Graphic interface NOT re-initialized")

# %% Initialize


# %% Import
data = pd.read_csv(matlab_csv_filepath, sep=';')
print(data.columns.values)

# %% Analyze
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

traces_len = data.trace_id.max() + 1
intersects = np.ones([traces_len, 1]) * np.nan
slopes = np.ones([traces_len, 1]) * np.nan

for trace_id in trange(1000 + 0 * traces_len):
    analyzed_count += 1
    cur_trace = data[data.trace_id == trace_id]

    frames_len = cur_trace.count().Frame

    # print(frames_len)
    # A trace must contain at least `expected_slope_length_frames` frames
    if frames_len < expected_slope_length_frames:
        skipped_too_short += 1
        short_trace_lengths.append(frames_len)
        continue

    start_frame = cur_trace.Frame.iloc[0]
    # Identify the first maximum fluo for the given trace
    max_frame = cur_trace.Frame[cur_trace.polymerases == cur_trace.polymerases.max()].iloc[0]

    # # If discontinuous traces, skip and log
    # if cur_trace.Frame.iloc[-1] - cur_trace.Frame.iloc[0] + 1 != cur_trace.count().Frame:
    #     skipped_discontinuous += 1
    #     continue

    # There must be at least `expected_slope_length_frames` before the maximum
    earlier_frames = cur_trace[cur_trace.Frame <= max_frame]
    if earlier_frames.count().Frame < expected_slope_length_frames:
        skipped_short_slope += 1
        short_slope_lengths.append(earlier_frames.count().Frame)
        continue

    # # The first frame must not be a maximum
    # if max_frame == start_frame:
    #     skipped_starts_with_the_maximum += 1
    #     continue

    # print(earlier_frames)

    # Get the possible locations for the start of the fixed-frame count fit
    earlier_frames = earlier_frames.iloc[0:(-expected_slope_length_frames + 1)]

    # Find the last minimum among these locations
    min_frame = earlier_frames[earlier_frames.polymerases ==
                               earlier_frames.polymerases.min()].Frame.iloc[-1]

    # Extract frames with the potential slope between the last minimum and the maximum
    potential_slope_frames = cur_trace[(cur_trace.Frame >= min_frame)
                                       & (cur_trace.Frame <= max_frame)]
    # slope_frames_len = potential_slope_frames.count().Frame
    # print([trace_id, min_frame, max_frame])
    # print(potential_slope_frames)

    # # If slope length is shorter than expected_slope_length_frames, skip and log
    # if slope_frames_len < expected_slope_length_frames:
    #     skipped_short_slope += 1
    #     continue

    slope, intersect, residual = identify_shift_and_slope(potential_slope_frames.time_min.values,
                                                          potential_slope_frames.polymerases.values, potential_slope_frames.nc.iloc[0])
    # print([slope, intersect, residual])

    # Drop if slope is negative
    if slope < 0:
        skipped_negative_slope += 1
        negative_slopes.append(slope)
        continue

    # Store
    slopes[trace_id] = slope
    intersects[trace_id] = intersect
    processed_successfully += 1

# Slope detection statistics
print('%i traces analyzed' % analyzed_count)
print('%i traces successfully processed' % processed_successfully)
print('%i traces skipped because too short. Median length: %.i' %
      (skipped_too_short, np.median(short_trace_lengths)))
# print('%i traces skipped because started with the maximum' % skipped_starts_with_the_maximum)
print('%i traces skipped because of short slope. Median slope length: %.i' %
      (skipped_short_slope, np.median(short_slope_lengths)))
print('%i traces skipped because negative slope. Median slope: %.2f' %
      (skipped_negative_slope, np.median(negative_slopes)))
# print('%i traces skipped because discontinuous' % skipped_discontinuous)


# %%
print(slopes.tolist())
print(intersects.tolist())
# print(data[data.trace_id == 56])
print(np.nanmedian(data.APpos))
