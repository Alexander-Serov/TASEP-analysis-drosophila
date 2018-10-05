import matplotlib.pyplot as plt
import os
from reinit_folder import reinit_folder
from tqdm import trange

xy_folder = r'.\xy_plots'


def plot_one_embryo(data, dataset_id):

    # Filter data
    cur_data = data[data.dataset_id == dataset_id].groupby(by='trace_id').first()
    dataset_name = cur_data.dataset.iloc[0]
    # print(cur_data)
    # print(cur_data)

    x_coord = 'AP'
    # x_coord = 'xPos'

    fig = plt.figure(11, clear=True)

    if x_coord == 'AP':
        plt.scatter(cur_data.MeanAP, cur_data.yPosMean)
        plt.xlabel('AP')
        plt.xlim([0, 1])
    elif x_coord == 'xPos':
        plt.scatter(cur_data.xPosMean, cur_data.yPosMean)
        plt.xlabel('xPos')

    plt.ylabel('yPosMean')
    str_title = '%s, id=%i' % (dataset_name, dataset_id)
    plt.title(str_title)
    # plt.show()

    figpath = f'xy_plot_{dataset_id}.png'
    figpath = os.path.join(xy_folder, figpath)
    plt.savefig(figpath)


def plot_xy_in_all_embryos(data):
    reinit_folder(xy_folder)

    dataset_id_max = data.dataset_id.max()

    for dataset_id in trange(dataset_id_max + 1):
        plot_one_embryo(data, dataset_id)
