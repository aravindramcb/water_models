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


def get_average_bottleneck(tunnels_def: list, simulation_results: str, tt_results: str, save_loc: str):
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
    save_name = os.path.join(save_loc, "average_bottleneck_main_figure.csv")
    avg_df.to_csv(save_name)
    return avg_df


def get_helix_distance(tunnel: str, simulation_results: str, save_loc: str) -> pd.DataFrame:
    distance_from_csv = pd.read_csv(f'/mnt/gpu/dean/md/{tunnel}_openings.csv', index_col='Unnamed: 0')
    directories = os.listdir(simulation_results)
    directories.sort(key=lambda x: (x.split("_")[0][:-1], x))
    distance_from_csv = distance_from_csv.reindex(columns=directories)
    models = ['opc', 'tip3p', 'tip4pew']
    groups = ["TCG0", "TCG1", "TCG2", "TCG3", "TCG4"]

    # get mean distance of group, regardless of models
    g = 0
    mean_df = pd.DataFrame()
    for i in range(5):
        column_name = groups[i]
        _df = distance_from_csv.iloc[:, g:g + 15]
        _mean_series = pd.Series([val for sublist in _df.values for val in sublist])
        mean_df[column_name] = _mean_series
        g += 15
    save_name = os.path.join(save_loc, "helix_distance_main_figure.csv")
    mean_df.to_csv(save_name)
    return mean_df


def load_csv_file(consolidated_csv_file):
    print("make sure you modify load_csv_file's usecols parameter to load correct columns")
    consolidated_df = pd.read_csv(consolidated_csv_file,usecols=[1,2,3,4,5])
    return consolidated_df


def save_defaultdict_to_csv(defaultdict, filename):
    import csv
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=defaultdict.keys())
        writer.writeheader()
        for key, value in defaultdict.items():
            writer.writerow({key: value})


def get_transit_time(tt_results, sim_results, group_definition, type: str = 'combined', save_loc: str = None):
    """
    Gets time spent by waters in the tunnel. Splits by models & groups
    :param type: type of event to get, it can be entry, release or combined
    :param tt_results: Transport Tools results location
    :param sim_results: Simulation results location
    :param group_definition: scid:[super_cluster numbers]
    :return: Dataframe per group
    """
    result = s7.get_transit_time(tt_results, sim_results, group_definition, debug=False)
    if type == 'combined':
        values = result[0]
    elif type == 'entry':
        values = result[1]
    elif type == 'release':
        values = result[2]
    else:
        print("Specify correct type - combined, entry, release")
    models = ['opc', 'tip3p', 'tip4pew']
    x_labels = ["TCG0", "TCG1", "TCG2", "TCG3", "TCG4"]
    names = ["1", "1.4", "1.8", "2.4", "3"]

    # Split into models + group
    i, j, k = (0, 5, 10)
    median_dict = defaultdict(pd.DataFrame)
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
        median_dict[x_labels[row]] = pd.DataFrame(_median_df).T
    save_name = os.path.join(save_loc, type + "_transit_time_median.csv")
    save_defaultdict_to_csv(median_dict, save_name)
    return median_dict


def save_to_obj(tuple, filename):
    import pickle
    with open(filename, 'wb') as f:
        pickle.dump(tuple, f)


def load_from_obj(filename):
    import pickle
    with open(filename, 'rb') as f:
        tuple = pickle.load(f)
    return tuple


def plot_waters_per_frame(tt_results: str, sim_results: str, tunnels_definition: dict, save_loc: str):
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
    saved_object = os.path.join(save_loc, "plot_waters_per_frames_fetched_frames.obj")
    if not os.path.isfile(saved_object):
        fetched_frames = t_events.get_transit_time(tt_results=tt_results, simulation_results=sim_results,
                                               groups_definitions=tunnels_definition,
                                               frame_numbers=True)
        # Save the fetched_frames for easy future plotting
        save_file_name = os.path.join(save_loc, "plot_waters_per_frames_fetched_frames.obj")
        save_to_obj(fetched_frames, save_file_name)

    fetched_frames = load_from_obj(saved_object)

    overall_color = sns.color_palette('deep', 3)
    sns.set(style="white", context="paper", font_scale=1)
    fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(8, 6), dpi=300)

    names = ["1.4", "1.8", "2.4", "3"]
    models = ["opc", "tip3p", "tip4pew"]
    events = ['Entry&Release', 'Entry', 'Release']
    for event_type in range(3):  # 0= entry&release, 1=entry, 2=release
        for group in range(4):  # 5 groups
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
                    average_water_per_frame = np.average(list(frames.values()))
                    plot_data.append(average_water_per_frame)

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

            ax[event_type, group].set_xticklabels([])
            ax[event_type, group].set_ylim(0, 6.5)
        # ax[event_type, 2].set_title(f"WATERS PER FRAME {events[event_type]} - P1 tunnel\n", fontweight='bold')
        ax[event_type, 0].set_ylabel("Avg. Num waters")
        ax[event_type, 0].set_xlabel("TCG1")
        ax[event_type, 1].set_xlabel("TCG2")
        ax[event_type, 2].set_xlabel("TCG3")
        ax[event_type, 3].set_xlabel("TCG4")

    plt.tight_layout()
    save_figure_name = os.path.join(save_loc, "water_per_frame.png")
    plt.savefig(save_figure_name)


def figure_two(caver_bottleneck, helix_distance, save_location: str = None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    groups = ["TCG0", "TCG1", "TCG2", "TCG3", "TCG4"]
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
    # hd.set_title("A)  Helix-Helix distance of P1 tunnel")
    hd.set_xticklabels(groups)
    hd.set_ylabel("Distance (Å)")
    hd.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Bottleneck radii
    ax1 = fig.add_subplot(ga[1, :])
    bp = sns.boxplot(data=caver_bottleneck, ax=ax1, boxprops=boxprops, whiskerprops=whiskerprops,
                     capprops=capprops, flierprops=flierprops)
    # bp.set_title("Average bottleneck radii of P1 tunnel")
    bp.set_xticklabels(groups)
    bp.set_ylabel("Bottleneck radii")

    save_name = os.path.join(save_location, "figure2.png")
    plt.tight_layout()
    plt.savefig(save_name)


def plot_water_transit_time(tunnels_definition, tt_results, simulation_results, save_location):
    sns.set(style="white", context="paper", font_scale=1)
    fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(8, 6), dpi=300)
    plt.suptitle("Transit time median")
    groups = ["TCG1", "TCG2", "TCG3", "TCG4"]
    # y_limits = [110, 95, 120, 95]
    y_limits =[120,120,120,120]
    obj = os.path.join(save_location,"plot_water_transit_time_rt_full.obj")
    if not os.path.exists(obj):
        rt_full = get_transit_time(tt_results, simulation_results, tunnels_definition, type='combined',
                                   save_loc=save_location)
        # Save the fetched_frames for easy future plotting
        save_file_name = os.path.join(save_location, "plot_water_transit_time_rt_full.obj")
        save_to_obj(rt_full, save_file_name)
    else:
        rt_full = load_from_obj(obj)

    # Entry and release
    for i in range(4):
        data = rt_full[groups[i]]
        data = data.apply(lambda x: x * 10, axis=0)
        print(groups[i], "Entry&Release\n", data)
        ts = sns.barplot(data=data, ax=ax[0, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts.set_xticklabels([])
        ts.set_xlabel(groups[i])
        # TO set uniform y axes per group
        ts.set_ylim(0, y_limits[i])
        i += 1

    # ax[0, 1].set_title("Water transit time medians [Entry + Release] - P1 Tunnel ".upper(), fontweight='bold',loc='left')
    ax[0, 0].set_ylabel("Time (ps)", fontweight="bold")

    # Entry
    rt_entry = get_transit_time(tt_results, simulation_results, tunnels_definition, type='entry',
                                save_loc=save_location)
    for i in range(4):
        data = rt_entry[groups[i]]
        data = data.apply(lambda x: x * 10, axis=0)
        print(groups[i], "Entry\n", data)
        ts_entry = sns.barplot(data=data, ax=ax[1, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts_entry.set_xticklabels([])
        ts_entry.set_xlabel(groups[i])
        # TO set uniform y axes per group
        ts_entry.set_ylim(0, y_limits[i])
        i += 1
    # ax[1, 1].set_title("Water transit time medians [Entry] - P1 Tunnel ".upper(), fontweight='bold',loc='left')
    ax[1, 0].set_ylabel("Time (ps)", fontweight="bold")

    # Release
    rt_release = get_transit_time(tt_results, simulation_results, tunnels_definition, type='release',
                                  save_loc=save_location)
    for i in range(4):
        data = rt_release[groups[i]]
        data = data.apply(lambda x: x * 10, axis=0)
        print(groups[i], "Release\n", data)
        ts_release = sns.barplot(data=data, ax=ax[2, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts_release.set_xticklabels([])
        ts_release.set_xlabel(groups[i])
        # TO set uniform y axes per group
        ts_release.set_ylim(0, y_limits[i])
        i += 1
    # ax[2, 1].set_title("Water transit time medians [Release] - P1 Tunnel".upper(), fontweight='bold',loc='center')
    ax[2, 0].set_ylabel("Time (ps)", fontweight="bold")

    save_file = os.path.join(save_location, "transit_time_median.png")
    plt.tight_layout(h_pad=2)
    plt.savefig(save_file)

    print(f"saved to {save_file}")
    plt.close()


def tt_events(tt_events, save_location, normailzed=None):
    sns.set(style="white", context="paper", font_scale=1.3)
    fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(10, 3), dpi=300)
    groups = [ "TCG1", "TCG2", "TCG3", "TCG4"]
    # Transport events
    j = 3
    tt_axes = []
    for i in range(4):
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

    tt_axes[2].set_title("Transport Events - P1\n".upper(), fontweight="bold",loc='left')
    if normailzed is not None:
        tt_axes[0].set_ylabel(f"Events - Normalized")
        save_name = os.path.join(save_location, "tt_events_norm.png")
    else:
        tt_axes[0].set_ylabel("Number of transported \n water molecules")
        save_name = os.path.join(save_location, "tt_events.png")

    plt.tight_layout()
    plt.savefig(save_name)


def plot_percent_event_occurrence(tt_results: str, sim_results: str, tunnels_def: dict, save_location: str):
    sns.set(style="white", context="paper", font_scale=0.8)
    fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(8, 6), dpi=300)
    groups = ["TCG1", "TCG2", "TCG3", "TCG4"]
    plt.suptitle("ratio of conductivity")
    events = s7.fraction_events_occurrence(tt_results, sim_results, tunnels_def)

    # save all loaded events to a file
    save_name= os.path.join(save_location,"plot_percent_event_occurrence_events.obj")
    save_to_obj(events,save_name)

    entry_release = events[1:5]
    for i in range(4):
        data = entry_release[i]
        data = data.apply(lambda x: x * 100, axis=0)
        # print(groups[i], "Entry&Release\n", data)
        ts = sns.barplot(data=data, ax=ax[0, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts.set_xticklabels([])
        ts.set_xlabel(groups[i])
        ts.set_ylim(0, 100)
        i += 1

    # insets because the first 2 plots have very small percentages

    data = entry_release[0]
    data = data.apply(lambda x: x * 100, axis=0)
    left, bottom, width, height = [0.25, 0.4, 0.5, 0.4]
    inset_axes = ax[0, 0].inset_axes([left, bottom, width, height])
    inset_plot = sns.barplot(data=data, ax=inset_axes, width=0.5, errorbar='se', capsize=.1, linewidth=1,
                             errwidth=1)
    inset_plot.set_xticks([])
    # inset_plot.set_yticks([])
    inset_plot.set_xlabel('')
    inset_plot.set_ylabel('')
    inset_plot.set_ylim(0, 4)

    # ax[0, 2].set_title("% of frames per simulation involving Entry + Release Events - P1 Tunnel ".upper(),
    #                    fontweight='bold')
    ax[0, 0].set_ylabel("% Frames", fontweight="bold")

    entry = events[6:10]
    for i in range(4):
        data = entry[i]
        data = data.apply(lambda x: x * 100, axis=0)
        print(groups[i], "Entry\n", data)
        ts_entry = sns.barplot(data=data, ax=ax[1, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts_entry.set_xticklabels([])
        ts_entry.set_xlabel(groups[i])
        ts_entry.set_ylim(0,100)
        i += 1
    # ax[1, 2].set_title("% of frames per simulation involving Entry Events - P1 Tunnel ".upper(), fontweight='bold')
    ax[1, 0].set_ylabel("% Frames", fontweight="bold")

    # inset
    data = entry[0]
    data = data.apply(lambda x: x * 100, axis=0)
    left, bottom, width, height = [0.25, 0.4, 0.5, 0.4]
    inset_axes = ax[1, 0].inset_axes([left, bottom, width, height])
    inset_plot = sns.barplot(data=data, ax=inset_axes, width=0.5, errorbar='se', capsize=.1, linewidth=1,
                             errwidth=1)
    inset_plot.set_xticks([])
    # inset_plot.set_yticks([])
    inset_plot.set_xlabel('')
    inset_plot.set_ylabel('')
    inset_plot.set_ylim(0, 4)

    release = events[11:15]
    for i in range(4):
        data = release[i]
        data = data.apply(lambda x: x * 100, axis=0)
        print(groups[i], "Release\n", data)
        ts_release = sns.barplot(data=data, ax=ax[2, i], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ts_release.set_xticklabels([])
        ts_release.set_xlabel(groups[i])
        ts_release.set_ylim(0,100)
        i += 1
    # ax[2, 2].set_title("% of frames per simulation involving Release Events- P1 Tunnel ".upper(), fontweight='bold')
    ax[2, 0].set_ylabel("% Frames", fontweight="bold")

    # inset
    data = release[0]
    data = data.apply(lambda x: x * 100, axis=0)
    left, bottom, width, height = [0.25, 0.4, 0.5, 0.4]
    inset_axes = ax[2, 0].inset_axes([left, bottom, width, height])
    inset_plot = sns.barplot(data=data, ax=inset_axes, width=0.5, errorbar='se', capsize=.1, linewidth=1,
                             errwidth=1)
    inset_plot.set_xticks([])
    # inset_plot.set_yticks([])
    inset_plot.set_xlabel('')
    inset_plot.set_ylabel('')
    inset_plot.set_ylim(0, 4)

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
    simulation_results = "/data/aravindramt/dean/md/simulations/"
    save_location = "/home/aravind/PhD_local/dean/figures/main_images/"
    consolidated_csv_file = '/home/aravind/PhD_local/dean/figures/transport_tools/p1_only.csv'
    bottleneck_csv ="/home/aravind/PhD_local/dean/figures/main_images/average_bottleneck_main_figure.csv"
    helix_csv = "/home/aravind/PhD_local/dean/figures/main_images/helix_distance_main_figure.csv"

    # bottleneck = get_average_bottleneck(tunnels_def=P1, simulation_results=simulation_results, tt_results=tt_results
    #                                     ,save_loc=save_location)
    # helix = get_helix_distance("p1", simulation_results,save_loc=save_location)
    # tt = load_csv_file(consolidated_csv_file)
    bottleneck = load_csv_file(bottleneck_csv)
    helix = load_csv_file(helix_csv)
    figure_two(bottleneck, helix, save_location)
    # tt_events(tt, save_location)
    # plot_water_retention_time(rt, save_location, normailzed="bygroup")

    # NEW !!
    # plot_water_transit_time(tunnels_definition=main_tunnel, tt_results=tt_results,
    #                         simulation_results=simulation_results, save_location=save_location)

    # plot_waters_per_frame(tt_results=tt_results, sim_results=simulation_results, tunnels_definition=main_tunnel,
    #                       save_loc=save_location)
    # plot_percent_event_occurrence(tt_results, simulation_results, main_tunnel, save_location)


if __name__ == '__main__':
    main()
