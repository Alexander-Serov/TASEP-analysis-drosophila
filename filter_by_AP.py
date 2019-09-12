import pandas as pd


def filter_by_AP(data_in):
    # data = data_in.copy()

    # data.groupby[by = 'dataset_id']
    datasets = set(data_in.dataset)
    data = pd.DataFrame()

    for dataset in datasets:
        dataset_data = data_in[data_in.dataset == dataset]
        gene = dataset_data.gene.iloc[0]
        max_AP = dataset_data.ap_mean.max()
        median_AP = dataset_data.ap_mean.median()

        if gene == 'hb':
            # do not include the transition region of width 0.05 at the right side
            AP_trans_region_width = 0.05
            filtered_data = dataset_data[(dataset_data.ap_mean <= max_AP - AP_trans_region_width)
                                         & (dataset_data.ap_mean >= 0.05)]
            # print(max_AP)
        elif gene == 'kn':
            AP_central_region_width = 0.04
            filtered_data = dataset_data[(dataset_data.ap_mean >= median_AP - AP_central_region_width / 2) &
                                         (dataset_data.ap_mean <= median_AP + AP_central_region_width / 2)]
        elif gene == 'sn':
            # No filtering for the moment. Update if necessary
            # AP_central_region_width = 0.04
            filtered_data = dataset_data
        data = pd.concat([data, filtered_data])

    return data
