

import pandas as pd
from tqdm import tqdm


def filter_by_AP(data_in):
    """
    This function contains the code used for filtering fluorescence traces by their AP positions

    For hb, we find the AP corresponding to maximum expression APmax, and then take all traces in the interval [0.05; APmax - 0.05].

    For kn, the expressing AP region moves with time. So instead of doing filtering over the whole trace altogether, we do frame by frame filtering. We do not need the connectivity of the traces, because they are anyway averaged out in a single trace for each data set.
    The kn expressing region is rather small, so the filtering works by finding the median AP (mAP) of all expressing nuclei for each individual frame and then keeping only traces withing mAP +- 0.02.

    For sn, no filtering is performed.
    """

    # Constants
    knirps_half_width_AP = 0.02
    hb_AP_margin = 0.05

    datasets = set(data_in.dataset)
    data = pd.DataFrame()

    for dataset in tqdm(datasets, desc='Parsing data sets'):
        # Take one data set
        dataset_data = data_in[data_in.dataset == dataset]
        gene = dataset_data.gene.iloc[0]

        if gene == 'hb':
            # Localize the maximum and keep the AP region described above
            max_AP = dataset_data.ap_mean.max()
            filtered_data = dataset_data[
                (dataset_data.ap_mean <= max_AP - hb_AP_margin) &
                (dataset_data.ap_mean >= hb_AP_margin)]
            data = pd.concat([data, filtered_data])

        elif gene == 'kn':
            frames = set(dataset_data.frame)
            # Individually process each frame
            for frame in frames:
                frame_data = dataset_data[dataset_data.frame == frame]
                median_AP = frame_data.ap.median()
                filtered_data = frame_data[
                    (frame_data.ap >= median_AP - knirps_half_width_AP) &
                    (frame_data.ap <= median_AP + knirps_half_width_AP)]
                data = pd.concat([data, filtered_data])

        elif gene == 'sn':
            filtered_data = dataset_data
            data = pd.concat([data, filtered_data])
    print('Data filtered!')

    return data
