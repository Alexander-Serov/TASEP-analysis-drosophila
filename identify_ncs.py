

import numpy as np
import os
import pandas as pd
from save_series_plot import save_series_plot
from tqdm import trange

from constants import mins_per_frame, min_nc13_length_minutes, max_nc13_length_minutes, nc13_folder, max_reasonable_polymerase_number


def identify_ncs(data):
    """
    Locate the start of nc13, nc14 and the end of nc13 to correctly measure the slopes
    """

    datasets_len = data.dataset_id.max() + 1
    ncs_locations = pd.DataFrame(
        columns=['dataset_id', 'nc13_start', 'nc13_end', 'nc14_start'], index=range(0, datasets_len))

    skipped_polymerases_too_high = 0

    # %% Identification of starts and ends of the ncs 13 and 14
    # dataset_id = 2
    for dataset_id in trange(datasets_len):
        # dataset_id = 29
        [dataset_id]

        start_nc13_frame = np.nan
        end_nc13_frame = np.nan
        start_nc14_frame = np.nan
        bl_search_nc13 = True

        # Select data
        dataset_nc_data = data[(data.dataset_id == dataset_id)]

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

        # If the dataset features polymerase numbers higher than 200, that must be a data error
        # The data will be rechecked
        if avg_dataset_nc_data.polymerases.max() >= max_reasonable_polymerase_number:
            [start_nc13_frame, end_nc13_frame,
                start_nc14_frame] = np.nan * np.asarray([1, 1, 1])
            skipped_polymerases_too_high += 1

            # Save ncs locations
        ncs_locations.loc[dataset_id] = [dataset_id,
                                         start_nc13_frame, end_nc13_frame, start_nc14_frame]

        # Plot
        fit_interval = [0, 0]  # np.asarray(fit_frames) * mins_per_frame * 0
        vbars = np.asarray([start_nc13_frame, end_nc13_frame, start_nc14_frame]) * mins_per_frame
        filename = 'slopes_dataset_%i.png' % (dataset_id)
        filepath = os.path.join(nc13_folder, filename)
        save_series_plot(x=avg_dataset_nc_data.time_min,
                         y=avg_dataset_nc_data.polymerases, a=0, b=0, fit_interval=fit_interval, vbars=vbars, filename=filepath, dataset_name=dataset_name, dataset_id=dataset_id)
        # return

    # Report errors
    if skipped_polymerases_too_high > 0:
        print("Warning: %i datasets skipped due to polymerase values >= %.0f" %
              (skipped_polymerases_too_high, max_reasonable_polymerase_number))

    return ncs_locations
