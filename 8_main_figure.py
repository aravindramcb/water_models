# -*- coding: utf-8 -*-
__author__ = 'Aravind Selvaram Thirunavukarasu'
__email__ = 'arathi@amu.edu.pl, aravind1233@gmail.com'

import os
from collections import defaultdict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from libs import transport_events_analysis as s7, time_evolution_bottleneck as s4
import seaborn as sns


def get_average_bottleneck(tunnels_def: list, simulation_results: str, tt_results: str):
    sc_details_loc = os.path.join(tt_results, "data", "super_clusters", "details", "initial_super_cluster_details.txt")
    original_ids_dict = s4.get_orig_caver_id(req_sc_ids=tunnels_def, initial_sc_details_txt=sc_details_loc,
                                             simulation_results_dir=simulation_results)
    bottlenecks = s4.get_bottleneck_radii(scid_orig_caver_ID=original_ids_dict,
                                          sim_results_location=simulation_results)
    avg_df = s4.process_bottleneck(bottlenecks)
    # print average and std of each group
    to_print_average = avg_df.mean(axis=0)
    to_print_std = avg_df.std(axis=0)
    print(f"Average: \n--------------\n{to_print_average}\n-------------\nStandard Deviation\n--------\n{to_print_std}")
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


def get_transit_time(tt_results, sim_results, group_definiton, type: str = 'combined'):
    """
    Gets time spent by waters in the tunnel. Splits by models & groups
    :param type: type of event to get, it can be entry, release or combined
    :param tt_results: Transport Tools results location
    :param sim_results: Simulation results location
    :param group_definiton: scid:[super_cluster numbers]
    :return: Dataframe per group
    """
    result = s7.get_transit_time(tt_results, sim_results, group_definiton, debug=False)
    if type == 'combined':
        values = result[0]
    elif type == 'entry':
        values = result[1]
    elif type == 'release':
        values = result[2]
    else:
        print("Specify correct type - combined, entry, release")
    models = ['opc', 'tip3p', 'tip4pew']
    x_labels = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]
    names = ["1", "1.4", "1.8", "2.4", "3"]

    # Split into models + group
    i, j, k = (0, 5, 10)
    median_df = defaultdict(pd.DataFrame)
    for row in range(5):  # 5 groups
        _median_df = []

        for col in range(3):  # 3 water models
            _median_simulations = []
            current_model = models[col]

            for sim_number in range(1, 6):  # 5 simulations - starting from 1-5
                sim_name = f"{names[row]}A_{current_model}_{sim_number}"
                _median = np.median(values[sim_name])
                _median_simulations.append(_median)
            _median_df.append(_median_simulations)
            i += 5
        median_df[x_labels[row]] = pd.DataFrame(_median_df).T
    return median_df


def plot_waters_per_frame(tt_results: str, sim_results: str, tunnels_definition: dict):
    """
    consolidated 1x5 plot of average and standard deviation of number of events per frame. The avg & std will
    imply how much waters are present per frame. This will say that in different groups on an average there are x amount
    of waters present in y'th frame.(old description, to be updated)
    :param tt_results: TT results location
    :param sim_results: Simulation results location
    :param tunnels_definition: Tunnels definition (dict)
    :return:
    """
    from libs import transport_events_analysis as t_events
    import numpy as np
    from collections import Counter
    fetched_frames = t_events.get_transit_time(tt_results=tt_results, simulation_results=sim_results,
                                               groups_definitions=tunnels_definition,
                                               frame_numbers=True)
    overall_color = sns.color_palette('deep', 3)
    sns.set(style="white", context="paper", font_scale=1)
    fig, ax = plt.subplots(nrows=3, ncols=5, figsize=(8, 6), dpi=300)

    names = ["1", "1.4", "1.8", "2.4", "3"]
    models = ["opc", "tip3p", "tip4pew"]
    events = ['Entry&Release', 'Entry', 'Release']
    for event_type in range(3):  # 0= entry&release, 1=entry, 2=release
        for group in range(5):  # 5 groups
            group_df = pd.DataFrame()
            print(f"{events[event_type]} Group {group + 1}")

            for model in models:
                group_name = f"{model}_{names[group]}"

                # Get frame numbers for 5 simulations in the current group
                data = []
                for sim in range(1, 6):
                    sim_name = f"{names[group]}A_{model}_{sim}"
                    frames_in_current_event_type = fetched_frames[event_type]
                    _frames = frames_in_current_event_type[sim_name]
                    data.append(_frames)
                # Get the average frames per simulation for every group
                plot_data = []  # This will have average number of frames per simulation
                for sim in data:
                    frames = Counter(sim)
                    average_event_per_frame = np.average(list(frames.values()))
                    plot_data.append(average_event_per_frame)

                # set color
                if "opc" in group_name:
                    color = overall_color[0]
                elif "tip3p" in group_name:
                    color = overall_color[1]
                else:
                    color = overall_color[2]

                # create dataframes to plot
                if "opc" in group_name:
                    df = pd.DataFrame(plot_data)
                    group_df[group_name] = df
                elif "tip3p" in group_name:
                    df = pd.DataFrame(plot_data)
                    group_df[group_name] = df
                else:
                    df = pd.DataFrame(plot_data)
                    group_df[group_name] = df

            print(group_df)

            # barplot
            sns.barplot(data=group_df, ax=ax[event_type, group], palette=overall_color, errorbar='sd',
                        capsize=.1, linewidth=1, errwidth=1)

            # Violin Plot if needed
            # sns.violinplot(data=group_df,ax=ax[group],palette=overall_color)
            # ax[group].get_children()[5].set_color('red')
            # for i in [2,6,10]:
            #     ax[group].get_children()[i].set_color('r')

            ax[event_type, group].set_xticklabels([])
            ax[event_type, group].set_ylim(0, 6.5)
        ax[event_type, 2].set_title(f"WATERS PER FRAME {events[event_type]} - P1 tunnel\n", fontweight='bold')
        ax[event_type, 0].set_ylabel("Avg. Num waters")
        ax[event_type, 0].set_xlabel("Group 1")
        ax[event_type, 1].set_xlabel("Group 2")
        ax[event_type, 2].set_xlabel("Group 3")
        ax[event_type, 3].set_xlabel("Group 4")
        ax[event_type, 4].set_xlabel("Group 5")
    plt.tight_layout()
    plt.savefig(f"/home/aravind/PhD_local/dean/figures/main_images/water_per_frame.png")


def figure_two(caver_bottleneck, helix_distance, save_location: str = None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    groups = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]
    sns.set(style='ticks')
    fig = plt.figure(figsize=(15, 6), dpi=300, constrained_layout=False)
    sns.set_context(context="paper", font_scale=2)
    boxprops = {'facecolor': 'silver', 'edgecolor': 'black', 'linewidth': 0.5}
    whiskerprops = {'linewidth': 0.5}
    capprops = {'linewidth': 0.5}
    flierprops = {"marker": ".", "markersize": 1}
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
                     capprops=capprops, flierprops=flierprops)
    bp.set_title("B) Average bottleneck radii of P1 tunnel")
    bp.set_xticklabels(groups)
    bp.set_ylabel("Bottleneck (Å)")

    save_name = os.path.join(save_location, "figure2.png")
    plt.tight_layout()
    plt.savefig(save_name)


def plot_water_transit_time(tunnels_definition, tt_results, simulation_results, save_location):
    sns.set(style="white", context="paper", font_scale=0.8)
    fig, ax = plt.subplots(nrows=3, ncols=5, figsize=(8, 6), dpi=300)
    groups = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]
    y_limits = [50,110,95,120,95]
    rt_full = get_transit_time(tt_results, simulation_results, tunnels_definition, type='combined')

   # Entry and release
    for i in range(5):
        data = rt_full[groups[i]]
        data = data.apply(lambda x: x * 10, axis=0)
        print(groups[i], "Entry&Release\n", data)

        ts = sns.barplot(data=data, ax=ax[0, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts.set_xticklabels([])
        ts.set_xlabel(groups[i])
        # TO set uniform y axes per group
        ts.set_ylim(0,y_limits[i])
        i += 1

    ax[0, 2].set_title("Water transit time medians [Entry + Release] - P1 Tunnel ".upper(), fontweight='bold')
    ax[0, 0].set_ylabel("Time (ps)", fontweight="bold")

    # Entry
    rt_entry = get_transit_time(tt_results, simulation_results, tunnels_definition, type='entry')
    for i in range(5):
        data = rt_entry[groups[i]]
        data = data.apply(lambda x: x * 10, axis=0)
        print(groups[i], "Entry\n", data)
        ts_entry = sns.barplot(data=data, ax=ax[1, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts_entry.set_xticklabels([])
        ts_entry.set_xlabel(groups[i])
        # TO set uniform y axes per group
        ts_entry.set_ylim(0, y_limits[i])
        i += 1
    ax[1, 2].set_title("Water transit time medians [Entry] - P1 Tunnel ".upper(), fontweight='bold')
    ax[1, 0].set_ylabel("Time (ps)", fontweight="bold")

    # Release
    rt_release = get_transit_time(tt_results, simulation_results, tunnels_definition, type='release')
    for i in range(5):
        data = rt_release[groups[i]]
        data = data.apply(lambda x: x * 10, axis=0)
        print(groups[i], "Release\n", data)
        ts_release = sns.barplot(data=data, ax=ax[2, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts_release.set_xticklabels([])
        ts_release.set_xlabel(groups[i])
        # TO set uniform y axes per group
        ts_release.set_ylim(0, y_limits[i])
        i += 1
    ax[2, 2].set_title("Water transit time medians [Release] - P1 Tunnel ".upper(), fontweight='bold')
    ax[2, 0].set_ylabel("Time (ps)", fontweight="bold")

    save_file = os.path.join(save_location, "transit_time_median.png")
    plt.tight_layout()
    plt.savefig(save_file)
    plt.close()


def tt_events(tt_events, save_location, normailzed=None):
    sns.set(style="white", context="paper", font_scale=1.3)
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(10, 3), dpi=300)
    groups = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]
    # Transport events
    j = 0
    tt_axes = []
    for i in range(5):
        _df = tt_events.iloc[:, j:j + 3]
        # print(_df)
        if normailzed == 'bygroup':
            max_value_bygroup = max(_df.max())
            axs = sns.barplot(data=_df.apply(lambda x: x / max_value_bygroup, axis=0), ax=ax[i], width=0.5,
                              errorbar='se',
                              capsize=.1, linewidth=1, errwidth=1)
            # axs.set_ylim(0, 1)
            print(_df.apply(lambda x: x / max_value_bygroup, axis=0))

        elif normailzed == 'whole':
            max_value_whole = max(tt_events[col].max() for col in tt_events.columns)
            axs = sns.barplot(data=_df.apply(lambda x: x / max_value_whole, axis=0), ax=ax[i], width=0.5, errorbar='se',
                              capsize=.1, linewidth=1, errwidth=1)
            axs.set_ylim(0, 1)
            print(_df.apply(lambda x: x / max_value_whole, axis=0))

        else:
            print(_df)
            axs = sns.barplot(data=_df, ax=ax[i], width=0.5, errorbar='se',
                              capsize=.1, linewidth=1, errwidth=1)

        axs.set_xticklabels([])
        axs.set_xlabel(groups[i])
        # axs.set_ylim(0, 1)
        tt_axes.append(axs)
        j += 3

    tt_axes[2].set_title("Transport Events - P1\n".upper(), fontweight="bold")
    if normailzed is not None:
        tt_axes[0].set_ylabel(f"Events - Normalized")
        save_name = os.path.join(save_location, "tt_events_norm.png")
    else:
        tt_axes[0].set_ylabel("Number of Events")
        save_name = os.path.join(save_location, "tt_events.png")

    plt.tight_layout()
    plt.savefig(save_name)


def plot_percent_event_occurrence(tt_results: str, sim_results: str, tunnels_def: dict, save_location: str):
    sns.set(style="white", context="paper", font_scale=0.8)
    fig, ax = plt.subplots(nrows=3, ncols=5, figsize=(8, 6), dpi=300)
    groups = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5"]
    events = s7.fraction_events_occurrence(tt_results, sim_results, tunnels_def)
    entry_release = events[0:5]
    for i in range(5):
        data = entry_release[i]
        data = data.apply(lambda x: x * 100, axis=0)
        # print(groups[i], "Entry&Release\n", data)
        ts = sns.barplot(data=data, ax=ax[0, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts.set_xticklabels([])
        ts.set_xlabel(groups[i])
        ts.set_ylim(0, 100)
        i += 1

    # insets becaue the first 2 plots have very small percentages
    # for i in range(2):
    #     data = entry_release[i]
    #     data = data.apply(lambda x: x * 100, axis=0)
    #     left, bottom, width, height = [0.25, 0.4, 0.5, 0.4]
    #     inset_axes = ax[0, i].inset_axes([left, bottom, width, height])
    #     inset_plot = sns.barplot(data=data, ax=inset_axes, width=0.5, errorbar='se', capsize=.1, linewidth=1,
    #                              errwidth=1)
    #     inset_plot.set_xticks([])
    #     # inset_plot.set_yticks([])
    #     inset_plot.set_xlabel('')
    #     inset_plot.set_ylabel('')
    #     inset_plot.set_ylim(0, 4)

    ax[0, 2].set_title("% of frames per simulation involving Entry + Release Events - P1 Tunnel ".upper(),
                       fontweight='bold')
    ax[0, 0].set_ylabel("% Frames", fontweight="bold")

    entry = events[5:10]
    for i in range(5):
        data = entry[i]
        data = data.apply(lambda x: x * 100, axis=0)
        print(groups[i], "Entry\n", data)
        ts_entry = sns.barplot(data=data, ax=ax[1, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts_entry.set_xticklabels([])
        ts_entry.set_xlabel(groups[i])
        i += 1
    ax[1, 2].set_title("% of frames per simulation involving Entry Events - P1 Tunnel ".upper(), fontweight='bold')
    ax[1, 0].set_ylabel("% Frames", fontweight="bold")

    release = events[10:15]
    for i in range(5):
        data = release[i]
        data = data.apply(lambda x: x * 100, axis=0)
        print(groups[i], "Release\n", data)
        ts_release = sns.barplot(data=data, ax=ax[2, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts_release.set_xticklabels([])
        ts_release.set_xlabel(groups[i])
        i += 1
    ax[2, 2].set_title("% of frames per simulation involving Release Events- P1 Tunnel ".upper(), fontweight='bold')
    ax[2, 0].set_ylabel("% Frames", fontweight="bold")

    # save_file = os.path.join(save_location, "fraction_frames_events_inset.png")
    save_file = os.path.join(save_location, "percent_frames_events.png")
    plt.tight_layout()
    plt.savefig(save_file)
    plt.close()


def main():
    # Set the variables
    # ------------------------------------------------------- #
    P1 = [1, 2, 5, 7, 12, 30, 31]
    main_tunnel = {"P1": [1, 2, 5, 7, 12, 30, 31]}
    P2 = [3, 4, 6, 8, 11, 16, 25, 27, 41, 43, 44, 50, 58]
    P3 = [10]
    tt_results = '/data/aravindramt/dean/tt/tt_0_9_5/'
    simulation_results = "/data/aravindramt/dean/tt/minimal_data"
    save_location = "/home/aravind/PhD_local/dean/figures/main_images/"
    consolidated_csv_file = '/home/aravind/PhD_local/dean/figures/transport_tools/p1_only.csv'

    bottleneck = get_average_bottleneck(tunnels_def=P1, simulation_results=simulation_results, tt_results=tt_results)
    # helix = get_helix_distance("p1", simulation_results)
    # tt = get_tt_events(consolidated_csv_file)
    # figure_two(bottleneck, helix, save_location)
    # tt_events(tt, save_location)
    # plot_water_retention_time(rt, save_location, normailzed="bygroup")

    # NEW !!
    # plot_water_transit_time(tunnels_definition=main_tunnel, tt_results=tt_results,
    #                         simulation_results=simulation_results, save_location=save_location)
    #
    # plot_waters_per_frame(tt_results=tt_results, sim_results=simulation_results, tunnels_definition=main_tunnel)
    # plot_percent_event_occurrence(tt_results, simulation_results, main_tunnel, save_location)


if __name__ == '__main__':
    main()
