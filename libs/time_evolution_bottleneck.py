# Python script to plot the time evolution of bottlenecks per simulation, grouped by model and group name.
# -*- coding: utf-8 -*-

# Author: Aravind Thirunavukarasu
# Email : arathi@amu.edu.pl, aravind1233@gmail.com
import json
import os
from collections import defaultdict

import pandas
import pandas as pd
from pandas import DataFrame


def get_orig_caver_id(req_sc_ids:list, initial_sc_details_txt:str, simulation_results_dir:str):
    """
    For the given SuperCluster IDs, this will get the highest priority corresponding original caver cluster ID from
    initial_super_cluster_details.txt
    :param initial_sc_details_txt: The file location of initial_super_clusters_details.txt
    :type req_sc_ids: SuperCluster IDs for your group, example in my case P1 is formed by the following SuperClusters
    ['1', "2", "5", "7", "12", "16", "30", "31"]
    """
    dirs = [d for d in os.listdir(simulation_results_dir) if os.path.isdir(os.path.join(simulation_results_dir, d))]
    result = {d: [] for d in dirs}
    with open(initial_sc_details_txt) as f:
        content = f.read()
        clusters = content.split(
            "--------------------------------------------------------------------------------------------"
            "----------------------------")
        for cluster in clusters:
            cluster = cluster.strip()
            lines = cluster.split("\n")
            if lines[0].startswith("Supercluster ID"):
                super_cluster_id = int(lines[0].split()[-1])
                if super_cluster_id in req_sc_ids:
                    for line in lines[6:]:
                        if line.startswith("from"):
                            sim_id, values = line.split(":")
                            values = [x.strip() for x in values.split(',') if x.strip()]
                            sim_id = sim_id.split()[1]
                            print(line)
                            try:
                                # convert values to int
                                values = [int(x) for x in values]
                                if sim_id not in result:
                                    result[sim_id] = values
                                else:
                                    result[sim_id].extend(values)
                            except ValueError:
                                pass
                else:
                    continue
        if any(not v for v in result.values()):
            print("Some simulations have missing corresponding caver clusters assigned in TransportTools")
            missing_values = []
            for key, value in result.items():
                if not value:
                    missing_values.append(key)
            print(missing_values, sep="|")
        else:
            print("All simulations have corresponding caver clusters assigned in TransportTools ")
        reduced_result = {k: min(v) if v else None for k, v in result.items()}
        return reduced_result


def get_bottleneck_radii(scid_orig_caver_ID: dict, sim_results_location: str):
    """
    Gets the bottleneck radii for given dict of simulation_ID:original_caver_cluster_id from tunnel_characteristics.csv from
    the original caver result for that simulation_ID
    :param scid_orig_caver_ID: Dictionary containing simulation_ID:original_caver_cluster_id
    :param sim_results_location: Simulation directory locaton
    :return: Dataframe of sim_ID as column name and bottleneck radii as values for the simulation frames.
    """
    bottleneck_radius = pd.DataFrame()
    # print("The assigned caver clusters are: \n")
    for sim_id, cav_cluster_id in scid_orig_caver_ID.items():
        # print(sim_id, "->", cav_cluster_id)
        csv_loc = 'caver_analyses/final_clustering/analysis/tunnel_characteristics.csv'
        csv_file = os.path.join(sim_results_location, sim_id, csv_loc)
        tc = pd.read_csv(csv_file, usecols=[1, 5])
        tc_filtered = tc.loc[(tc[" Tunnel cluster"] == cav_cluster_id)]
        bottleneck_radius[sim_id] = (tc_filtered[" Bottleneck radius"].reset_index(drop=True))
    bl_sorted_df = (bottleneck_radius.reindex(sorted(bottleneck_radius.columns), axis=1))
    bl_sorted_df.to_csv("/home/aravind/PhD_local/dean/figures/bottlenecks/time_evolution/bottleneck_per_sim.csv")
    return bl_sorted_df


def plot_bottlenecks(bottleneck_dataframe: DataFrame, group_name: str):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    bl_sorted_df = bottleneck_dataframe

    # # Prepare data
    col_names = bl_sorted_df.columns.values.tolist()
    # df_melt = pd.melt(bl_sorted_df.reset_index(), id_vars=['index'], value_vars=col_names)
    # df_melt.columns = ['index', 'Sim_ID', 'Bottleneck']
    #
    # # Group data
    # avg_1 = pd.melt(bl_sorted_df.iloc[:, 30:45].reset_index(), id_vars='index', value_vars=col_names[30:45])
    # avg_1.columns = ['index', 'Sim_ID', 'Whole']
    # avg_1_4 = pd.melt(bl_sorted_df.iloc[:, 0:15].reset_index(), id_vars='index', value_vars=col_names[0:15])
    # avg_1_4.columns = ['index', 'Sim_ID', 'Whole']
    # avg_1_8 = pd.melt(bl_sorted_df.iloc[:, 15:30].reset_index(), id_vars='index', value_vars=col_names[15:30])
    # avg_1_8.columns = ['index', 'Sim_ID', 'Whole']
    # avg_2_4 = pd.melt(bl_sorted_df.iloc[:, 45:60].reset_index(), id_vars='index', value_vars=col_names[45:60])
    # avg_2_4.columns = ['index', 'Sim_ID', 'Whole']
    # avg_3 = pd.melt(bl_sorted_df.iloc[:, 60:75].reset_index(), id_vars='index', value_vars=col_names[60:75])
    # avg_3.columns = ['index', 'Sim_ID', 'Whole']
    #
    # avg_df = pd.concat([avg_1.loc[:, 'Whole'], avg_1_4.loc[:, 'Whole'], avg_1_8.loc[:, 'Whole'],
    #                     avg_2_4.loc[:, 'Whole'], avg_3.loc[:, 'Whole']], axis=1)
    # avg_df.columns = ['1A', '1.4A', '1.8A', '2.4A', '3A']
    #
    # # Concat data per group
    # df_1 = pd.concat([avg_1.loc[:, 'Whole'], bl_sorted_df.iloc[:, 30:45]], axis=1)
    # df_1_4 = pd.concat([avg_1_4.loc[:, 'Whole'], bl_sorted_df.iloc[:, 0:15]], axis=1)
    # df_1_8 = pd.concat([avg_1_8.loc[:, 'Whole'], bl_sorted_df.iloc[:, 15:30]], axis=1)
    # df_2_4 = pd.concat([avg_2_4.loc[:, 'Whole'], bl_sorted_df.iloc[:, 45:60]], axis=1)
    # df_3 = pd.concat([avg_3.loc[:, 'Whole'], bl_sorted_df.iloc[:, 60:75]], axis=1)
    #
    # df without whole

    df_1 = bl_sorted_df.iloc[:, 30:45]
    df_1_4 = bl_sorted_df.iloc[:, 0:15]
    df_1_8 = bl_sorted_df.iloc[:, 15:30]
    df_2_4 = bl_sorted_df.iloc[:, 45:60]
    df_3 = bl_sorted_df.iloc[:, 60:75]

    x_labels = [col_names[i].split(sep='_', maxsplit=1)[1] for i in range(15)]
    fig, axes = plt.subplots(3, 2, figsize=(30, 25), dpi=150)
    plt.subplots_adjust(hspace=0.3, top=1)
    sns.set()
    plt.suptitle(f"BOTTLENECK RADII OF TUNNEL {group_name}", fontsize=20, fontweight='bold')
    # plt.suptitle("BOTTLENECK RADII OF P3 TUNNEL", fontsize=20, fontweight='bold',y=0.98)
    box1 = sns.boxplot(data=df_1, ax=axes[0, 0], color='b')
    box1.set_xlabel("TCG0", fontsize=20, fontweight='bold')
    # box1.artists[0].set_facecolor('grey')

    box2 = sns.boxplot(data=df_1_4, ax=axes[0, 1], color='g')
    box2.set_xlabel("TCG1", fontsize=20, fontweight='bold')
    # box2.artists[0].set_facecolor('grey')

    box3 = sns.boxplot(data=df_1_8, ax=axes[1, 0], color='r')
    box3.set_xlabel("TCG2", fontsize=20, fontweight='bold')
    box3.set_ylabel("Bottleneck radii (Å)", fontsize=20, fontweight='bold')
    # box3.artists[0].set_facecolor('grey')

    box4 = sns.boxplot(data=df_2_4, ax=axes[1, 1], color='c')
    box4.set_xlabel("TCG3", fontsize=20, fontweight='bold')
    # box4.artists[0].set_facecolor('grey')

    box5 = sns.boxplot(data=df_3, ax=axes[2, 0], color='m')
    box5.set_xlabel("TCG4", fontsize=20, fontweight='bold')
    # box5.artists[0].set_facecolor('grey')
    x_tick_location = np.arange(15)
    for ax in axes.flatten():
        ax.set_ylim(bottom=0.8, top=4)
        ax.set_xticks(x_tick_location)
        ax.set_xticklabels(x_labels, rotation=30, fontsize=20)
        ax.tick_params(axis='y', labelsize=20)
        sns.set_theme(style='darkgrid')
    color_pal = {'1A': 'b', '1.4A': 'g', '1.8A': 'r', '2.4A': 'c', '3A': 'm'}
    # box6 = sns.boxplot(data=avg_df, palette=color_pal)
    # box6.set_xlabel("Overall", fontsize=20, fontweight='bold')
    # box6.set_xticklabels(box6.get_xticklabels(),["TCG0","TCG1","TCG2","TCG3","TCG4"])
    plt.tight_layout(pad=1.8)
    plt.savefig(f"/home/aravind/PhD_local/dean/figures/bottlenecks/time_evolution/{group_name}_manuscript.png")


def process_bottleneck(bottleneck_dataframe):
    bl_sorted_df = bottleneck_dataframe

    # Prepare data
    col_names = bl_sorted_df.columns.values.tolist()
    df_melt = pd.melt(bl_sorted_df.reset_index(), id_vars=['index'], value_vars=col_names)
    df_melt.columns = ['index', 'Sim_ID', 'Bottleneck']

    # Group data
    avg_1 = pd.melt(bl_sorted_df.iloc[:, 30:45].reset_index(), id_vars='index', value_vars=col_names[30:45])
    avg_1.columns = ['index', 'Sim_ID', 'Whole']
    avg_1_4 = pd.melt(bl_sorted_df.iloc[:, 0:15].reset_index(), id_vars='index', value_vars=col_names[0:15])
    avg_1_4.columns = ['index', 'Sim_ID', 'Whole']
    avg_1_8 = pd.melt(bl_sorted_df.iloc[:, 15:30].reset_index(), id_vars='index', value_vars=col_names[15:30])
    avg_1_8.columns = ['index', 'Sim_ID', 'Whole']
    avg_2_4 = pd.melt(bl_sorted_df.iloc[:, 45:60].reset_index(), id_vars='index', value_vars=col_names[45:60])
    avg_2_4.columns = ['index', 'Sim_ID', 'Whole']
    avg_3 = pd.melt(bl_sorted_df.iloc[:, 60:75].reset_index(), id_vars='index', value_vars=col_names[60:75])
    avg_3.columns = ['index', 'Sim_ID', 'Whole']

    avg_df = pd.concat([avg_1.loc[:, 'Whole'], avg_1_4.loc[:, 'Whole'], avg_1_8.loc[:, 'Whole'],
                        avg_2_4.loc[:, 'Whole'], avg_3.loc[:, 'Whole']], axis=1)
    avg_df.columns = ['1A', '1.4A', '1.8A', '2.4A', '3A']

    return avg_df


def plot_bottlenecks_overview(tunnels_def: list, group_names: list, simulation_results: str, sc_details_loc: str,
                              save_loc: str):
    import seaborn as sns
    import matplotlib.pyplot as plt

    bottlenecks_of_groups = defaultdict(pandas.DataFrame)
    i = 0
    for group in group_names:
        original_ids_dict = get_orig_caver_id(req_sc_ids=tunnels_def[i], initial_sc_details_txt=sc_details_loc,
                                              simulation_results_dir=simulation_results)
        bottlenecks = get_bottleneck_radii(scid_orig_caver_ID=original_ids_dict,
                                           sim_results_location=simulation_results)
        bottlenecks_of_groups[group] = bottlenecks
        i += 1

    average_df = defaultdict(pandas.DataFrame)
    for group in group_names:
        avg_df = process_bottleneck(bottlenecks_of_groups[group])
        average_df[group] = avg_df
    sns.set(style='whitegrid')
    fig, axes = plt.subplots(nrows=1, ncols=3, dpi=300, figsize=(10, 4))
    color_pal = {'1A': 'b', '1.4A': 'g', '1.8A': 'r', '2.4A': 'c', '3A': 'm'}
    box1 = sns.boxplot(ax=axes[0], data=average_df['P1'], palette=color_pal, linewidth=0.5, fliersize=0.5)
    box1.set_xlabel("P1", fontweight="bold")
    box1.set_ylabel("Bottleneck radii (Å)")
    box1.set_ylim(0.7, 3.7)
    box2 = sns.boxplot(ax=axes[1], data=average_df['P2'], palette=color_pal, linewidth=0.5, fliersize=0.5)
    box2.set_xlabel("P2", fontweight="bold")
    box2.set_ylim(0.7, 3.7)
    box3 = sns.boxplot(ax=axes[2], data=average_df['P3'], palette=color_pal, linewidth=0.5, fliersize=0.5)
    box3.set_xlabel("P3", fontweight="bold")
    box3.set_ylim(0.7, 3.7)
    plt.suptitle("TIME EVOLUTION OF BOTTLENECKS", fontweight="bold")
    save_location = os.path.join(save_loc + "overall.png")
    plt.tight_layout()
    plt.savefig(save_location)


if __name__ == '__main__':
    # Set the variables
    P1 = [1, 2, 5, 7, 12, 30, 31]
    P2 = [3, 4, 6, 8, 11, 16, 25, 27, 41, 43, 44, 50, 58]
    P3 = [10]
    sc_details_file_loc = "/data/aravindramt/dean/tt/tt_0_9_5/data/super_clusters/details" \
                          "/initial_super_cluster_details.txt"
    simulation_results = "/data/aravindramt/dean/md/simulations/"
    save_location = "/home/aravind/PhD_local/dean/figures/bottlenecks/time_evolution/"

    # Get the original caver IDs for the given group name (P1,P2,P3) for all simulations
    original_ids_dict = get_orig_caver_id(req_sc_ids=P1, initial_sc_details_txt=sc_details_file_loc,
                                          simulation_results_dir=simulation_results)
    #
    # Get bottlenecks for the original Ids for all simulations
    bottlenecks = get_bottleneck_radii(scid_orig_caver_ID=original_ids_dict, sim_results_location=simulation_results)

    # Plot per group
    bottlenecks = pd.read_csv("/home/aravind/PhD_local/dean/figures/bottlenecks/time_evolution/bottleneck_per_sim.csv")
    plot_bottlenecks(bottleneck_dataframe=bottlenecks, group_name="P1")

    # Plot overall
    # plot_bottlenecks_overview(tunnels_def=[P1,P2,P3],group_names=["P1","P2","P3"],simulation_results=simulation_results,
    #                           sc_details_loc = sc_details_file_loc,save_loc=save_location)
