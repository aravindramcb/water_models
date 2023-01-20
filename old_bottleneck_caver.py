"""Script to check if the average bottlenecks of the tunnels are matching
with the bottlenecks of the snapshots"""
import numpy as np
import os
from matplotlib import pyplot as plt


# Get bottleneck from individual caver calculations
def parse_bottleneck_into_clusters(work_dir):
    """
    get the average,standard deviation and number of tunnels for a single epoch
    :param work_dir:
    :return:[bottleneck average, bottleneck std, number of tunnels]
    """
    bottleneck_avg = []
    bottleneck_std = []
    number_of_tunnels = []
    dir_list = ["1", "2", "3", "4", "5"]
    csv_dir = "/caver_analyses/final_clustering/analysis/"

    def parse_bottleneck(file_to_be_processed):
        """
        Gets the bottleneck radii and cluster number from tunnel_characterestics.csv file and filters only cluster 1
        :param file_to_be_processed:
        :return: average and standard deviation of the bottleneck radius of cluster 1
        """
        bottleneck = np.loadtxt(file_to_be_processed, skiprows=1, usecols=(5,), delimiter=',')
        cluster = np.loadtxt(file_to_be_processed, dtype=int, skiprows=1, usecols=(1,), delimiter=',')
        number_of_tunnels.append(bottleneck.size)
        filter_cluster = np.stack((bottleneck, cluster), axis=1)
        first_cluster_only = filter_cluster[filter_cluster[:, 1] < 2]  # Filter only cluster number <2
        return np.average(first_cluster_only[:, 0]), np.std(first_cluster_only[:, 0])

    # Iterate into directory and fetch csv
    for directory in dir_list:
        to_change = os.path.join(work_dir, (directory + csv_dir))
        final_clustering_dir = os.path.join(work_dir, (directory + '/caver_analyses/final_clustering'))
        try:
            os.chdir(final_clustering_dir)
        except NotADirectoryError:
            pass
        except FileNotFoundError:
            pass
        final_clustering_contents = os.listdir()

        if "analysis" in final_clustering_contents:
            analysis = os.chdir(to_change)

            if "tunnel_characteristics.csv" in os.listdir(analysis):
                # print(f'currently working on {directory}')
                avg_btlneck = parse_bottleneck('tunnel_characteristics.csv')
                bottleneck_std.append(avg_btlneck[1])
                bottleneck_avg.append(avg_btlneck[0])
        else:
            print(f'warning ! no analysis found in {directory}')
            pass
    return bottleneck_avg, bottleneck_std, number_of_tunnels


colour_list = [["blue", "orange", "green", "red", "purple"],
               ["purple", "red", "lawngreen", "orange", "blue"],
               ["aqua", "fuchsia", "dodgerblue", "pink", "coral"]]
models = ['opc', 'tip3p', 'tip4pew']
epoch_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
sim_list = [1,2,3,4,5]
expected = [1, 1.4, 1.8, 2.4, 3]


def plot_bottleneck_all_models(avg, std):
    rows, cols = 5, 5
    fig, ax = plt.subplots(5, 5, sharey='col', figsize=(25, 25))
    for system in [0, 1, 2]:
        for sim in range(cols):
            for epoch in range(rows):
                average = avg[system]  # Get average of a single system
                standard_deviation = std[system]  # Get std of the current system
                ax[sim, epoch].set_title(f"{epoch_list[epoch]}-{sim_list[sim]}")
                # lab = "Avg and std for"+ models[system]
                x_axes = sim_list[sim]
                y_axes = avg[system][epoch][sim]
                error = std[system][epoch][sim]
                ax[sim, epoch].errorbar(x_axes, y_axes, yerr=error, linestyle='None', marker="^", label=lab,
                                        color=colour_list[system][epoch])
                ax[sim, epoch].legend()
                # ax[sim,epoch].plot(x_axes)
        ax.set_label()
        # plt.show()
    fig.text(0.5, 0.025, 'Simulation Number', ha='center', weight='bold', fontsize=15, )
    fig.text(0.03, 0.5, 'Bottleneck radius ($\AA$)', weight='bold', fontsize=15, rotation='vertical')
    fig.suptitle("Bottleneck radius of caver cluster 1 for different water models", weight='bold', fontsize=20)
    plt.tight_layout(pad=9, w_pad=0, h_pad=4)
    plt.savefig(fname="/home/aravind/test2.png")


def plot_btlneck_combined(avg, std):
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(1, 5,figsize=(18, 8), sharey='col')
    for system in [0, 1, 2]:
        for sim in range(5):
            for epoch in range(5):
                # Get std of the current system
                custom_lines = [Line2D([0], [0], color=colour_list[0][epoch], marker='^', lw=0),
                                Line2D([0], [0], color=colour_list[1][epoch], marker='^', lw=0),
                                Line2D([0], [0], color=colour_list[2][epoch], marker='^', lw=0)]
                ax[epoch].set_title(f"{epoch_list[epoch]}")
                lab = "Avg and std for" + models[system]
                x_axes = sim_list[sim]
                y_axes = avg[system][epoch][sim]
                error = std[system][epoch][sim]
                ref = expected[sim]
                ax[epoch].set_ylim(bottom=0.5, top=3.5)

                ax[epoch].errorbar(x_axes, y_axes, yerr=error, linestyle='None', marker="^", label=lab,
                                   color=colour_list[system][epoch])
                ax[epoch].legend(custom_lines, ['OPC', 'TIP3P', 'TIP4PEW'])


    ax[0].plot(sim_list, [1,1,1,1,1])
    ax[1].plot(sim_list,[1.4,1.4,1.4,1.4,1.4])
    ax[2].plot(sim_list,[1.8,1.8,1.8,1.8,1.8])
    ax[3].plot(sim_list,[2.4,2.4,2.4,2.4,2.4])
    ax[4].plot(sim_list,[3,3,3,3,3])
                # ax[sim,epoch].plot(x_axes)
        # plt.show()
    fig.text(0.5, 0.025, 'Simulation Number', ha='center', weight='bold', fontsize=14)
    fig.text(0.01, 0.4, 'Bottleneck Radius ($\AA$)', weight='bold', fontsize=14, rotation='vertical')
    fig.suptitle("Average & Standard Deviation of different Bottlenecks for main tunnel", weight='bold', fontsize=16)
    # plt.title("Bottleneck radius of caver cluster 1 for different water models")
    plt.tight_layout(pad=3, w_pad=0, h_pad=2)
    plt.savefig(fname="/home/aravind/PhD_local/dean/distance_main_tunnel/combined.png")
    # plt.show()
    # plt.savefig(fname="/home/aravind/test2.png")
    # # print(bottleneck_std, bottleneck_avg)
    # print(number_of_tunnels)
    # # Plot graph with average and standard deviation
    # x = dir_list
    # y = bottleneck_avg
    # e = bottleneck_std
    # plt.title(subdir + " " + water_model)
    # expected = ((subdir[:-1] + ',') * 5).split(',', maxsplit=5)
    # plt.plot(x, list(map(float, expected[:-1])), label='Expected')
    # plt.errorbar(x, y, yerr=e, linestyle='None', marker='^', label='Average & Std')
    # plt.yticks(np.arange(0, 3.2, step=0.2))
    # plt.legend()

    # plt.close()


# top=0.856,
# bottom=0.136,
# left=0.076,
# right=0.974,
# hspace=0.14,
# wspace=0.231

# main

water_models = ['opc', 'tip3p', 'tip4pew']
avg_model, std_model, num_model = [], [], []
for model in water_models:
    root_dir = '/mnt/NAS/dean_water_models/md/' + model
    save_dir = root_dir + "/images/"
    subdir_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
    avg_epoch, std_epoch, numb_epoch = [], [], []
    for subdir in subdir_list:
        results_dir = '/mnt/NAS/dean_water_models/md/' + model + '/' + subdir
        avg, std, numb = parse_bottleneck_into_clusters(results_dir)
        avg_epoch.append([*avg])
        std_epoch.append([*std])
        numb_epoch.append([*numb])
    avg_model.append([*avg_epoch])
    std_model.append([*std_epoch])
    num_model.append([*numb_epoch])

plot_btlneck_combined(avg_model, std_model)
