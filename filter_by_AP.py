"""
For hb, the filtering is performed simply by the location along the AP position, the max. expression region.

For sn, no filtering for the moment.

For kn, the problem is that the expressing AP region moves with time. So instead of doing simple filtering, we re-evalute the expression median for each frame indpendently, and take the region +- 0.02 AP around it. I don't need to maintain connectivity of the traces since they are averaged out anyway later.
"""


import pandas as pd
from tqdm import tqdm

# Constants
knirps_half_width_AP = 0.02


def filter_by_AP(data_in):
    # data = data_in.copy()

    # data.groupby[by = 'dataset_id']
    datasets = set(data_in.dataset)
    data = pd.DataFrame()

    for dataset in tqdm(datasets, desc='Parsing data sets'):
        dataset_data = data_in[data_in.dataset == dataset]
        gene = dataset_data.gene.iloc[0]
        # median_AP = dataset_data.ap_mean.median()

        if gene == 'hb':
            # do not include the transition region of width 0.05 at the right side
            max_AP = dataset_data.ap_mean.max()
            AP_trans_region_width = 0.05
            filtered_data = dataset_data[(dataset_data.ap_mean <= max_AP - AP_trans_region_width)
                                         & (dataset_data.ap_mean >= 0.05)]
            data = pd.concat([data, filtered_data])

        elif gene == 'kn':
            frames = set(dataset_data.frame)
            # print(median_AP)
            for frame in frames:  # tqdm(frames, desc='Analyzing frames'):
                # AP_central_region_width = 0.04
                frame_data = dataset_data[dataset_data.frame == frame]
                median_AP = frame_data.ap.median()
                filtered_data = frame_data[
                    (frame_data.ap >= median_AP - knirps_half_width_AP) &
                    (frame_data.ap <= median_AP + knirps_half_width_AP)]
                data = pd.concat([data, filtered_data])

        elif gene == 'sn':
            # No filtering for the moment. Update if necessary
            # AP_central_region_width = 0.04
            filtered_data = dataset_data
            data = pd.concat([data, filtered_data])
    print('Data filtered!')

    return data
