

import numpy as np
import pandas as pd


def detect_y_regions_in_snail(data):
    """
    Detect the regions of maximal expression for the snail genes, where the averaging will be applied
    """

    y_bin_width = 10

    # Filter the snail gene
    data = data[data.gene_id == 2]

    datasets_ids = data.dataset_id.drop_duplicates()

    dataset_id = datasets_ids.iloc[0]

    dataset = data[data.dataset_id == dataset_id]
    y_min = dataset.yPos.min()
    y_max = dataset.yPos.max()

    n_bins = np.ceil(y_max / y_bin_width)
    y_bins = np.arange(0, y_bin_width * (n_bins + 1), y_bin_width)
    times = dataset.time_min.drop_duplicates().sort_values().values
    # print(times)

    # For each time frame, calculate an average per bin for this data set
    def get_bin_average(data, bin):
        y_bin_center = y_bins[bin - 1] + y_bin_width / 2
        bin_data = data[(data.yPos >= y_bins[bin - 1]) & (data.yPos <= y_bins[bin])]
        # print(bin_data)
        if bin_data.count().iloc[0] == 0:
            return None

        avg = bin_data.groupby(by="Frame").mean()[['time_min', 'polymerases']]
        # print(avg)
        # avg = avg.set_index('time_min')
        avg = avg.rename(columns={'polymerases': bin})
        return(avg)

    binned_data = pd.DataFrame(columns=times,
                               index=range(int(n_bins + 1)))
    # binned_data.index.name = 'time_min'
    # print(binned_data)
    bin = 10
    bin_avg = get_bin_average(dataset, bin)
    print(bin_avg)
    binned_data.iloc[bin] = bin_avg
    print(binned_data)

    return(binned_data)
