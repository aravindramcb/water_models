# -*- coding: utf-8 -*-
# python program to process bottleneck residues from overall results and comparative analysis results
# -*- coding: utf-8 -*-

__author__ = 'Aravind Selvaram Thirunavukarasu'
__email__ = 'arathi@amu.edu.pl, aravind1233@gmail.com'

import os
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def process_bottleneck(results_file:str, scids:list, frequency_cutoff: float = 0.2) -> dict:
    """
    Gets the bottleneck residues and frequency for the given SCIDs for the given results_file
    :param results_file: Results file containing bottleneck residues statistics from TransportTools results
    :param scids: SuperCluster IDs to parse bottlenecks
    :param frequency_cutoff: Cutoff to consider the residue number if it has to be processed.
    """
    # bottlenecks = defaultdict(dict)
    bottlenecks = []
    with open(results_file, "r") as infile:
        readfile = infile.readlines()
        for line in readfile[18:]:
            contents = line.split()
            contents = [item.strip(",") for item in contents]
            sc_id = int(contents[0])
            total_frames = contents[1]
            if total_frames == '-':
                continue
            else:
                res_freq = {entry.split(":")[0]: float(entry.split(":")[1]) for entry in contents[2:]}
                filtered_freq = {residue: frequency for residue, frequency in res_freq.items() if
                             frequency >= frequency_cutoff}
            if sc_id in scids:
                # bottlenecks[sc_id] = filtered_freq
                bottlenecks=bottlenecks+list(filtered_freq.keys())
    return bottlenecks


def get_comparative_results(comparative_resutls_location:str, scid:list):
    comparative_results = defaultdict(dict)
    if os.getcwd() != comparative_resutls_location:
        os.chdir(comparative_resutls_location)
    folders = [d for d in os.listdir() if os.path.isdir(d)]
    for folder in folders:
        _fname = "2-filtered_tunnels_statistics_bottleneck_residues.txt"
        _path = os.path.join(folder,_fname)
        if os.path.isfile(_path):
            btlnk = process_bottleneck(_path, scid)
            comparative_results[folder] = btlnk
    return comparative_results


def plot(scids: list, bottleneck:dict, save_location:str, group_name:str=None ):

    # sns.set_context(context="paper",rc={'figure.dpi':300})
    if len(scids) == 1:
        df = pd.DataFrame.from_dict(bottleneck, orient='index')
        ax_single = sns.barplot(data=df,color='black')
        ax_single.set_xlabel("Residue Number")
        ax_single.set_ylabel("Frequency")
        plt.title(label=f"Bottleneck residues in 0.2% or more for Tunnel-{scids[0]}")
        to_save = os.path.join(save_location,f"bottleneck_{scids}.png")
        plt.savefig(to_save)
        plt.close()
    else:
        for scid in scids:
            data = bottleneck[scid]
            df = pd.DataFrame.from_dict(data, orient='index').T
            ax_single = sns.barplot(data=df, edgecolor='black',hatch='//',color='white')
            ax_single.set_xlabel("Residue Number")
            ax_single.set_ylabel("Frequency")
            ax_single.set_ylim(0,1)
            plt.title(label=f"Bottleneck residues in 0.2% or more for Tunnel-{scid}")
            to_save = os.path.join(save_location, f"bottleneck_{scid}.png")
            if os.path.exists(save_location):
                plt.savefig(to_save)
            else:
                os.makedirs(save_location)
                plt.savefig(to_save)
            plt.close()

    # make full plot
    sns.set(style='white', rc={'figure.figsize':(11,8)})
    full_data = defaultdict(list)
    for scid in scids:
        for residue, frequency in bottleneck[scid].items():
            full_data[residue].append(frequency)
    average_data ={res:sum(freq)/len(freq) for res,freq in full_data.items()}
    data = pd.DataFrame.from_dict(average_data,orient='index')
    sorted_df = data.sort_values(0, ascending=False).T
    avg_ax = sns.barplot(data=sorted_df,edgecolor='black',color='b')
    avg_ax.set_xlabel("Residue Number")
    avg_ax.set_ylabel("Frequency")
    avg_ax.set_xticklabels(avg_ax.get_xticklabels(),rotation=90,size=10)
    if group_name is None:
        plt.title(f"Average Bottlenek Residues for the Group")
    else:
        plt.title(f"Average Bottlenek Residues for Tunnel {group_name}")
    avg_ax.set_ylim(0,1)
    to_save = os.path.join(save_location,f"overall.png")
    plt.savefig(to_save,dpi=300)
    # plt.show()
    plt.tight_layout()
    plt.close()

    # top ten
    sns.set(style='white', rc={'figure.figsize': (6, 5)})
    top_ten = sorted_df.iloc[:, 0:10]
    top_ten_ax = sns.barplot(data=top_ten, edgecolor='black', color='grey')
    top_ten_ax.set_xlabel("Residue Number")
    top_ten_ax.set_ylabel("Frequency")
    # top_ten_ax.set_xticklabels(top_ten_ax.get_xticklabels(), rotation=90, size=10)
    if group_name is None:
        plt.title(f"Top Ten Bottlenek Residues - Overall")
    else:
        plt.title(f"Top Ten Bottlenek Residues for {group_name}")
    top_ten_ax.set_ylim(0, 1)

    to_save = os.path.join(save_location, f"top_ten.png")
    plt.tight_layout()
    plt.savefig(to_save, dpi=300)
    plt.close()


def plot_top_ten_subplots(comparative_result:dict,group_name:str,save_loc:str) -> None:
    """
    Plot subplots of comparative analysis results, only for top 5 for the group.
    :param comparative_result: Dictionary contating comparative analysis results
    :param group_name: Name of the group being processed. eg, P1, P2, P3 ...
    :return: Plots figure
    """
    sns.set(style='white')
    fig,axes = plt.subplots(nrows=3,ncols=5,figsize =(11.69,8.27))
    fig.suptitle(f"Bottleneck residues for {group_name} - Top 5".upper())
    axes = axes.ravel()
    folder_list = sorted(list(comparative_result.keys()))
    for i,folder in enumerate(folder_list):
        full_data = defaultdict(list)
        current_result = comparative_result[folder]
        for scid in current_result:
            for residue, frequency in current_result[scid].items():
                full_data[residue].append(frequency)
        average_data = {res: sum(freq) / len(freq) for res, freq in full_data.items()}
        data = pd.DataFrame.from_dict(average_data,orient='index')
        sorted_df = data.sort_values(0, ascending=False).T
        top_ten = sorted_df.iloc[:, 0:5]
        ax = axes[i]
        color_list=sns.color_palette("tab10")
        sns.barplot(data=top_ten,ax=ax, edgecolor='black', color=color_list[2])
        # ax.set_xticklabels(ax.get_xticklabels(),fontsize = 10,rotation=45)
        ax.yaxis.set_tick_params(length=1)
        ax.set_ylim(0,1)
        ax.set_title(folder.upper())
    axes[12].set_xlabel("Residue Number".upper())
    axes[5].set_ylabel("Frequency of occurrence".upper())
    plt.tight_layout(w_pad=0.5)
    plt.savefig(save_loc+f"overall_{group_name}.png")


def plot_comparative(comparative_result:dict,save_location:str) -> None:
    if os.getcwd() is not save_location:
        os.chdir(save_location)
    if not os.path.exists("comparative_results"):
        os.mkdir("comparative_results")
    dir_list = list(comparative_result.keys())
    for folder in dir_list:
        full_path = os.path.join("comparative_results",folder)
        if not os.path.exists(full_path):
            os.mkdir(full_path)
        scids = list(comparative_result[folder].keys())
        plot(scids,comparative_result[folder],full_path,folder)


if __name__ == '__main__':

    result_file = "/data/aravindramt/dean/tt/tt_0_9_5_bottleneck/statistics/" \
                  "2-filtered_tunnels_statistics_bottleneck_residues.txt"
    p1 = [1, 2, 5, 7, 12, 30, 31]
    p2 = [3, 4, 6, 8, 11, 16, 25, 27, 41, 43, 44, 50, 58]
    p3 = [10]

    # From overall result - by tunnel
    save_location = "/home/aravind/PhD_local/dean/figures/bottlenecks/residues/p3"
    # df = process_bottleneck(result_file, p3)
    # plot(p3, df,save_location,group_name='P3')

    # Comparative result - by tunnel
    # save_location_comparative = "/home/aravind/PhD_local/dean/figures/bottlenecks/"
    comparative_analysis_results = "/data/aravindramt/result = {defaultdict: 15} defaultdict(<class 'dict'>, {'tip3p_2.4': ['146', '173', '138', '169', '242', '269', '270', '142', '165', '149', '172', '170', '103', '243', '268', '145', '141', '206', '168', '39', '38', '173', '142', '172', '145', '169', '141', '242', '149', '146', '138'... Viewdean/tt/tt_0_9_5_bottleneck/statistics/comparative_analysis"
    result = get_comparative_results(comparative_analysis_results,p1)
    print(result)
    # plot_comparative(result,save_location_comparative)

    # Subplots by tunnel
    # save_location_subplots = "/home/aravind/PhD_local/dean/figures/bottlenecks/"
    # plot_top_ten_subplots(result,group_name="P3",save_loc=save_location_subplots)

    # Subplots by model

