# -*- coding: utf-8 -*-
__author__ = 'Aravind Selvaram Thirunavukarasu'
__email__ = 'arathi@amu.edu.pl, aravind1233@gmail.com'

import os
from collections import defaultdict
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from libs import transport_events_analysis as s7, time_evolution_bottleneck as s4
import seaborn as sns

def get_average_bottleneck(tunnels_def: list, simulation_results: str, tt_results: str):
    sc_details_loc=os.path.join(tt_results,"data","super_clusters","details","initial_super_cluster_details.txt")
    original_ids_dict = s4.get_orig_caver_id(req_sc_ids=tunnels_def, initial_sc_details_txt=sc_details_loc,
                                             simulation_results_dir=simulation_results)
    bottlenecks = s4.get_bottleneck_radii(scid_orig_caver_ID=original_ids_dict,
                                          sim_results_location=simulation_results)
    avg_df = s4.process_bottleneck(bottlenecks)
    return avg_df


def get_helix_distance(tunnel: str, simulation_results) -> pd.DataFrame:
    distance_from_csv = pd.read_csv(f'/mnt/gpu/dean/md/{tunnel}_openings.csv', index_col='Unnamed: 0')
    directories = os.listdir(simulation_results)
    directories.sort(key=lambda x: (x.split("_")[0][:-1], x))
    distance_from_csv = distance_from_csv.reindex(columns=directories)
    models = ['opc', 'tip3p', 'tip4pew']
    groups = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]

    # get mean distance of group, regardless of models
    g = 0
    mean_df = pd.DataFrame()
    for i in range(5):
        column_name = groups[i]
        _df = distance_from_csv.iloc[:, g:g + 15]
        _mean_series = pd.Series([val for sublist in _df.values for val in sublist])
        mean_df[column_name] = _mean_series
        g += 15
    return mean_df


def get_tt_events(consolidated_csv_file):
    consolidated_df = pd.read_csv(consolidated_csv_file)
    return consolidated_df


def get_retention_time_overview(tt_results, sim_results, group_definiton):
    """
    Gets time spent by waters in the tunnel. Splits by models & groups
    :param tt_results: Transport Tools results location
    :param sim_results: Simulation results location
    :param group_definiton: scid:[super_cluster numbers]
    :return: Dataframe per group
    """
    result = s7.get_transit_time_single(tt_results, sim_results, group_definiton, debug=False)
    combined_dict = result[0]
    models = ['opc', 'tip3p', 'tip4pew']
    x_labels = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]
    names = ["1", "1.4", "1.8", "2.4", "3"]
    # Split into models + group
    i, j, k = (0, 5, 10)
    mean_df= defaultdict(pd.DataFrame)
    for row in range(5):
        _mean_df = []
        for col in range(3):
            current_model = models[col]
            column_name = f"{current_model}_{names[row]}"
            _mean_df.append(combined_dict[column_name])
            i += 5
        mean_df[x_labels[row]] = pd.DataFrame(_mean_df).T
    return mean_df




def figure_two(caver_bottleneck, helix_distance, save_location: str = None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    groups = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]
    sns.set(style='ticks')
    fig = plt.figure(figsize=(5, 3), dpi=300, constrained_layout=False)
    sns.set_context(context="paper",font_scale=0.8)
    boxprops = {'facecolor': 'silver', 'edgecolor': 'black', 'linewidth': 0.5}
    whiskerprops = {'linewidth': 0.5}
    capprops = {'linewidth': 0.5}
    flierprops = {"marker": ".","markersize":1}
    ga = fig.add_gridspec(nrows=2, ncols=3)

    # Helix-Helix distance
    ax2 = fig.add_subplot(ga[0, :])
    hd = sns.boxplot(data=helix_distance, ax=ax2, boxprops=boxprops, whiskerprops=whiskerprops, capprops=capprops,
                     flierprops=flierprops)
    hd.set_title("A)  Helix-Helix distance of P1 tunnel")
    hd.set_xticklabels(groups)
    hd.set_ylabel("Distance (Å)")
    hd.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Bottleneck radii
    ax1 = fig.add_subplot(ga[1, :])
    bp = sns.boxplot(data=caver_bottleneck, ax=ax1, boxprops=boxprops, whiskerprops=whiskerprops,
                     capprops=capprops,flierprops=flierprops)
    bp.set_title("B)  Average bottleneck radii of P1 tunnel")
    bp.set_xticklabels(groups)
    bp.set_ylabel("Bottleneck (Å)")

    save_name = os.path.join(save_location, "figure2.png")
    plt.tight_layout()
    plt.savefig(save_name)

def plot_water_retention_time(time_spent, save_location, normailzed=None):
    sns.set(style="white", context="paper", font_scale=0.8)
    fig, ax = plt.subplots(nrows=1, ncols=5,figsize=(8,2), dpi=300)
    groups = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]

    # Water retention
    for i in range(5):
        data = time_spent[groups[i]]
        print(data)
        if normailzed == 'whole':
            # Normalize by max value in all groups
            max_value = max(val.max().max() for val in time_spent.values())
            ts = sns.barplot(data=data.apply(lambda x: (x / max_value) * 20, axis=0), ax=ax[i], width=0.5, errorbar='se',
                             capsize=.1, linewidth=1, errwidth=1)
            print(data.apply(lambda x: (x / max_value)))
            # ts.set_ylim(0, 1)
        elif normailzed == 'bygroup':
            # Normalize by max value in group
            # max_value = max(data.max())
            # print(f"Max value = {max_value}")
            ts = sns.barplot(data=data.apply(lambda x: x / x.max(), axis=0), ax=ax[i], width=0.5, errorbar='se',
                             capsize=.1, linewidth=1, errwidth=1)
            print(data.apply(lambda x: x / x.max()))
            # ts.set_ylim(0, 1)
        else:
            ts = sns.barplot(data=data.apply(lambda x: x * 20,axis=0), ax=ax[i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
            # ts = sns.violinplot(data=data.apply(lambda x: x * 20,axis=0),ax=ax[i],width=0.5,linewidth=0.5)
            print(data.apply(lambda x: x * 20,axis=0))
        ts.set_xticklabels([])
        ts.set_xlabel(groups[i])
        i += 1
    # ax[0].set_ylim(0,20)
    ax[2].set_title("Water transit time - P1".upper(), fontweight='bold')
    if normailzed is not None:
        ax[0].set_ylabel(f"Time - Normalized {normailzed}", fontweight="bold")
        save_file = os.path.join(save_location, "retention_time_normalized.png")
    else:
        ax[0].set_ylabel("Time (ps)", fontweight="bold")
        save_file = os.path.join(save_location,"retention_time.png")
    plt.tight_layout()
    plt.savefig(save_file)


def tt_events(tt_events, save_location, normailzed=None):
    sns.set(style="white", context="paper", font_scale=0.8)
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(8,2),dpi=300)
    groups = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]
    # Transport events
    j = 0
    tt_axes = []
    for i in range(5):
        _df = tt_events.iloc[:, j:j + 3]
        # print(_df)
        if normailzed == 'bygroup':
            max_value_bygroup = max(_df.max())
            axs = sns.barplot(data=_df.apply(lambda x: x/max_value_bygroup,axis=0), ax=ax[i],width=0.5, errorbar='se',
                              capsize=.1, linewidth=1, errwidth=1)
            # axs.set_ylim(0, 1)
            print(_df.apply(lambda x: x/max_value_bygroup,axis=0))

        elif normailzed == 'whole':
            max_value_whole = max(tt_events[col].max() for col in tt_events.columns)
            axs = sns.barplot(data=_df.apply(lambda x: x / max_value_whole, axis=0), ax=ax[i], width=0.5, errorbar='se',
                              capsize=.1, linewidth=1, errwidth=1)
            axs.set_ylim(0, 1)
            print(_df.apply(lambda x: x/max_value_whole,axis=0))

        else:

            axs = sns.barplot(data=_df, ax=ax[i], width=0.5, errorbar='sd',
                              capsize=.1, linewidth=1, errwidth=1)

        axs.set_xticklabels([], fontsize=8)
        axs.set_xlabel(groups[i], fontsize=8)

        tt_axes.append(axs)
        j += 3

    tt_axes[2].set_title("Transport Events - P1".upper(),fontweight="bold")
    if normailzed is not None:
        tt_axes[0].set_ylabel(f"Events - Normalized {normailzed}", fontweight="bold")
        save_name = os.path.join(save_location, "tt_events_norm.png")
    else:
        tt_axes[0].set_ylabel("Number of Events", fontweight="bold")
        save_name = os.path.join(save_location, "tt_events.png")

    plt.tight_layout()
    plt.savefig(save_name)


def main():
    # Set the variables
    # ------------------------------------------------------- #
    P1 = [1, 2, 5, 7, 12, 30, 31]
    main_tunnel ={"P1": [1, 2, 5, 7, 12, 30, 31]}
    P2= [3, 4, 6, 8, 11, 16, 25, 27, 41, 43, 44, 50, 58]
    P3 = [10]
    tt_results = '/data/aravindramt/dean/tt/tt_0_9_5/'
    simulation_results = "/data/aravindramt/dean/tt/minimal_data"
    save_location = "/home/aravind/PhD_local/dean/figures/main_images/"
    consolidated_csv_file = '/home/aravind/PhD_local/dean/figures/transport_tools/p1_only.csv'

    bottleneck = get_average_bottleneck(tunnels_def=P1, simulation_results=simulation_results,tt_results=tt_results)
    helix = get_helix_distance("p1", simulation_results)
    rt = get_retention_time_overview(tt_results, simulation_results, {"P1": [1, 2, 5, 7, 12, 30, 31]})
    # tt = get_tt_events(consolidated_csv_file)

    figure_two(bottleneck, helix, save_location)
    # tt_events(tt,save_location,normailzed='whole')
    plot_water_retention_time(rt, save_location, normailzed="bygroup")


if __name__ == '__main__':
    main()
