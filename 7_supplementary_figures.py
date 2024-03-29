# -*- coding: utf-8 -*-
# Supplementary figures for haloalkane dehalogenase
__author__ = 'Aravind Selvaram Thirunavukarasu'
__email__ = 'arathi@amu.edu.pl, aravind1233@gmail.com'

import math

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
from libs import transport_events_analysis as tt_events
from libs import time_evolution_bottleneck as btlnk


def plot_consolidated_result(consolidated_csv_file: str, unassigned_csv: str, plot_normalized: bool = False):
    from matplotlib.patches import Patch
    colors = sns.color_palette('deep', 3)
    consolidated_df = pd.read_csv(consolidated_csv_file)
    i, j, k = (0, 20, 40)  #
    rows, cols = (5, 4)
    sns.set(style='white')
    fig, ax = plt.subplots(nrows=5, ncols=5, figsize=(11.69, 8.27), dpi=300, sharex='col')
    if plot_normalized is True:
        plt.suptitle("Events Scalled by TIP3P & Consolidated by tunnels and models", fontsize=15, fontweight='bold')
    else:
        plt.suptitle("Events Consolidated by tunnels and models", fontsize=15, fontweight='bold')
    for row in range(rows):
        for col in range(cols):
            _df = consolidated_df.iloc[:, [i, j, k]]
            if plot_normalized:
                # Normalize by TIP3P
                _df_sum = _df.sum()
                scaled_df = _df_sum / _df_sum.max()
                print(scaled_df)
                scaled_df.plot.bar(ax=ax[row, col], color=colors, sharey='row')
            else:
                # Consolidated with standard deviation
                print(_df)
                sns.barplot(data=_df, ax=ax[row, col], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
            ax[row, col].set_xticklabels([])
            i += 1
            j += 1
            k += 1
    a, b, c = (0, 5, 10)
    unassigned_df = pd.read_csv(unassigned_csv)
    # Plot unassigned
    for unas_row in range(rows):
        _unas_df = unassigned_df.iloc[:, [a, b, c]]
        print(_unas_df)
        sns.barplot(data=_unas_df, ax=ax[unas_row, 4], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
        ax[unas_row, 4].set_xticklabels([])
        a += 1
        b += 1
        c += 1
    plt.tight_layout(pad=3, w_pad=0.1, h_pad=0.1)
    opc_patch = Patch(color=colors[0], label='OPC')
    tip3p_patch = Patch(color=colors[1], label='TIP3P')
    tip4pew_patch = Patch(color=colors[2], label='TIP4P-Ew')
    ax[0, 4].legend([opc_patch, tip3p_patch, tip4pew_patch], ["OPC", "TIP3P", "TIP4P-Ew"],
                    bbox_to_anchor=(-0.9, 1.15, 0, 0), loc="lower right",
                    borderaxespad=1, ncol=3)
    ax[0, 0].set_title("P1", fontsize=10, fontweight='bold')
    ax[0, 1].set_title("P2", fontsize=10, fontweight='bold')
    ax[0, 2].set_title("P3", fontsize=10, fontweight='bold')
    ax[0, 0].set_ylabel("Group 1", fontsize=10)
    ax[1, 0].set_ylabel("Group 2", fontsize=10)
    ax[2, 0].set_ylabel("NUMBER OF FRAMES \n\n  Group 3", fontsize=10)
    ax[3, 0].set_ylabel("Group 4", fontsize=10)
    ax[4, 0].set_ylabel("Group 5", fontsize=10)
    ax[0, 3].set_title("Others", fontsize=10, fontweight='bold')
    ax[0, 4].set_title("Unassigned", fontsize=10, fontweight='bold')
    # fig.text(x=0.01, y=0.5, s="NUMBER OF SNAPSHOTS", rotation=90)
    fig.text(x=0.45, y=0.01, s="TUNNELS & MODELS")

    if plot_normalized:
        plt.savefig('/home/aravind/PhD_local/dean/figures/transport_tools/consolidated_tt_events_scaled.png')
    else:
        plt.savefig('/home/aravind/PhD_local/dean/figures/transport_tools/consolidated_tt_mean_se_events.png')
    # plt.show()


def plot_results_per_tunnel(tt_results: str, groups_definitions: dict, model: str, unassigned_csv: str):
    global sc_ids, index
    sns.set(style='white', font_scale=0.8)
    colors = {"opc": ['xkcd:dusk', 'xkcd:cool blue'], "tip3p": ['xkcd:dirt', 'xkcd:browny orange'],
              "tip4pew": ['xkcd:moss green', 'xkcd:seaweed']}
    fig, ax = plt.subplots(nrows=5, ncols=5, figsize=(11.69, 8.27), dpi=300, sharex='col')
    plt.suptitle(f"{model}".upper() + " Events", fontsize=15, fontweight='bold')
    rows, cols = (5, 4)
    i, j = (0, 5)
    x_label = [1, 2, 3, 4, 5]
    y_limits = [300,250,3700,12000,20000]
    comparative_analysis_results = os.path.join(tt_results, 'statistics/comparative_analysis')
    groups_scs = tt_events.get_scids_of_groups(comparative_analysis_results,
                                               groups_definitions)  # get corresponding SC per group
    tunnel_id = list(groups_definitions.keys()) + ["others"]
    opc_contents = list(groups_scs.keys())[0:5]
    tip3p_contents = list(groups_scs.keys())[5:10]
    tip4pew_contents = list(groups_scs.keys())[10:15]
    for row in range(rows):
        for col in range(cols):
            if model == "opc":
                group_name = opc_contents[row]
                sc_ids = groups_scs[group_name][0][tunnel_id[col]]
            if model == "tip3p":
                group_name = tip3p_contents[row]
                sc_ids = groups_scs[group_name][0][tunnel_id[col]]
            if model == "tip4pew":
                group_name = tip4pew_contents[row]
                sc_ids = groups_scs[group_name][0][tunnel_id[col]]
            entry_df, release_df = tt_events.get_entry_release_events(sc_ids, tt_results, model)
            data = entry_df[i:j]  # To plot entry per tunnel
            release_data = release_df[i:j]
            _combined_df = pd.concat([data, release_data], axis=1)
            _combined_df.reset_index(inplace=True, drop=True)
            _combined_df["Total"] = _combined_df.sum(axis=1)

            # using seaborn to plot
            sns.barplot(ax=ax[row, col], data=_combined_df, x=_combined_df.index, y='Total', color=colors[model][0])
            sns.barplot(ax=ax[row, col], data=_combined_df, x=_combined_df.index, y='Release_events',
                        color=colors[model][1])

            # # Modify x and y ticks
            ax[row, col].set_ylim(0,y_limits[row])
            ax[row, col].set_ylabel("")
            ax[row, col].set_xlabel("")
            ax[row, col].set_xticks([0, 1, 2, 3, 4])
            ax[row, col].set_xticklabels(x_label, rotation=0)

        i += 5
        j += 5

    # Giving explicit index related to Haloalkane Dehalogenase, modify for other proteins.
    if model == 'opc':
        index = 0
    elif model == 'tip3p':
        index = 25
    elif model == 'tip4pew':
        index = 50

    # Load unassigned data
    unassigned_df = pd.read_csv(unassigned_csv, index_col=0)

    # Plot unassigned
    for unas_col in range(5):
        _unas_df = unassigned_df.iloc[:, index:index + 5].T
        _unas_df["Total"] = _unas_df.sum(axis=1)
        sns.barplot(ax=ax[unas_col, 4], data=_unas_df, x=_unas_df.index, y='Total', color=colors[model][0])
        sns.barplot(ax=ax[unas_col, 4], data=_unas_df, x=_unas_df.index, y='Release', color=colors[model][1])
        ax[unas_col, 4].set_ylabel("")
        ax[unas_col, 4].set_xlabel("")
        ax[unas_col, 4].set_xticks([0, 1, 2, 3, 4])
        ax[unas_col, 4].set_xticklabels(x_label)
        index += 5

    # Generate Legend
    plt.tight_layout(pad=1.9, w_pad=0.1, h_pad=0.5)
    topbar = plt.Rectangle((0, 0), 1, 1, fc=colors[model][1], edgecolor='none')
    bottombar = plt.Rectangle((0, 0), 1, 1, fc=colors[model][0], edgecolor='none')
    ax[0, 4].legend([bottombar, topbar], ["ENTRY", "RELEASE"], bbox_to_anchor=(1.1, 1.15, 0, 0), loc="lower right",
                    borderaxespad=1, ncol=2)

    # Add more details
    ax[0, 0].set_title("P1", fontsize=15, fontweight='bold')
    ax[0, 1].set_title("P2", fontsize=15, fontweight='bold')
    ax[0, 2].set_title("P3", fontsize=15, fontweight='bold')
    ax[0, 3].set_title("Others", fontsize=15, fontweight='bold')
    ax[0, 4].set_title("Unassigned", fontsize=15, fontweight='bold')
    ax[0, 0].set_ylabel("TCG0", fontsize=15)
    ax[1, 0].set_ylabel("TCG1", fontsize=15)
    ax[2, 0].set_ylabel("TCG2", fontsize=15)
    ax[3, 0].set_ylabel("TCG3", fontsize=15)
    ax[4, 0].set_ylabel("TCG4", fontsize=15)
    ax[4, 2].set_xlabel("Simulation #", fontsize=15)
    plt.savefig(f'/home/aravind/PhD_local/dean/figures/transport_tools/combined_events_{model}.png')
    # plt.show()

def plot_transit_time_overview(tt_results, scids, sim_results, model:str= "overall", save_location:str=None):
    from collections import defaultdict
    sns.set(style='white', font_scale=0.8,context='paper')
    overall_color = sns.color_palette('deep', 3)
    colors = {"opc": ['xkcd:dusk', 'xkcd:cool blue'], "tip3p": ['xkcd:dirt', 'xkcd:browny orange'],
              "tip4pew": ['xkcd:moss green', 'xkcd:seaweed'],"overall":overall_color}

    times = defaultdict(pd.DataFrame)
    for scid in scids:
        result = tt_events.get_transit_time(tt_results=tt_results, scids=scids[scid], simulation_results=sim_results)
        rt = result[0]  # overall result
        times[scid] = rt
    tunnel_names = list(scids.keys())
    cols,rows = (3,5)
    fig,ax = plt.subplots(nrows=5, ncols=3, dpi=300)
    for col in range(cols):
        i, j, k = (0, 5, 10)
        for row in range(rows):
            tunnel_name = tunnel_names[col]
            entry=times[tunnel_name][0]
            release = times[tunnel_name][1]
            _combined = pd.concat([entry, release], axis=1, keys=["entry", "release"])
            # individulal sim / models
            opc = _combined.iloc[i:i+5,:]

            tip3p = _combined.iloc[j:j+5,:]
            tip4pew = _combined.iloc[k:k+5,:]
            # Combined inlet + release (mean)
            _combined_overall = _combined.mean(axis=1)
            overall_df_per_group = pd.concat([_combined_overall[i:i+5].reset_index(drop=True),
                                              _combined_overall[j:j+5].reset_index(drop=True),
                                              _combined_overall[k:k+5].reset_index(drop=True)],
                                             axis=1,keys=["OPC","TIP3P","TIP4P-Ew"])
            print(tunnel_name,model)
            if model == "opc":
                data = opc.reset_index(drop=True)
            elif model == "tip3p":
                data = tip3p.reset_index(drop=True)
            elif model == "tip4pew":
                data = tip4pew.reset_index(drop=True)
            else:
                data = overall_df_per_group
                plot = sns.barplot(data=data,ax=ax[row,col],width=0.5, errorbar='se', capsize=.1,
                                   linewidth=1, errwidth=1)
            print(data)
            if model != "overall":
                plot = data.plot(kind="bar", stacked=True,ax=ax[row,col],legend=False,color=colors[model])
            for axes in ax.flatten():
                axes.set_xticklabels([])
            i += 5
            j += 5
            k += 5
    plt.suptitle(f"Average retention frames per event - {model.upper()}",fontweight="bold")
    # annotate
    ax[0,0].set_title("P1 TUNNEL",fontweight='bold')
    ax[0, 1].set_title("P2 TUNNEL",fontweight='bold')
    ax[0, 2].set_title("P3 TUNNEL",fontweight='bold')
    for i in range(5):
        ax[i,0].set_ylabel(f"GROUP {i+1}",fontweight="bold")
    ax[2,0].set_ylabel(f"NUMBER OF FRAMES\n\nGROUP 3")
    # legend
    top_bar = plt.Rectangle((0, 0), 1, 1, fc=colors[model][1], edgecolor='none')
    bottombar = plt.Rectangle((0, 0), 1, 1, fc=colors[model][0], edgecolor='none')
    extra_bar = plt.Rectangle((0, 0), 1, 1, fc=colors["overall"][2], edgecolor='none')
    if model != "overall":
        ax[0, 1].legend([bottombar, top_bar], ["ENTRY", "RELEASE"], bbox_to_anchor=(1.1, 1.15, 0, 0), loc="lower right",
                        borderaxespad=1, ncol=2)
    else:
        ax[0,1].legend([bottombar,top_bar,extra_bar],["OPC","TIP3P","TIP4P-Ew"],bbox_to_anchor=(1.3, 1.15, 0, 0), loc="lower right",
                       borderaxespad=1, ncol=3)
    plt.tight_layout()
    if model == "overall":
        save_file = os.path.join(save_location,"overall.png")
    else:
        name = f"{model}.png"
        save_file = os.path.join(save_location,name)
    plt.savefig(save_file)
def histogram_count_transit_time(tt_results:str, sim_results:str, gd:dict):
    """
    Process data for histogram and plot transit time per event.
    1 frame = 10ps in my system
    transit time per event = number of frames per event * 10 for a single entry or release event
    :param tt_results: Location of Transport tools results
    :param sim_results: Location of simulation results
    :return: saves a plot of histogram in the specified location
    """
    frames = tt_events.get_transit_time(tt_results=tt_results, simulation_results=sim_results,
                                        groups_definitions=gd,
                                        frame_numbers=False)

    combined = frames[0]
    models = ['opc', 'tip3p', 'tip4pew']
    names = ["1.4", "1.8", "2.4", "3"]
    group_number = [2,3,4,5]
    overall_color = sns.color_palette('deep', 3)
    row, col = (4, 3)  # 5 groups and 3 models
    sns.set_context(context="paper", font_scale=1)
    fig,ax = plt.subplots(nrows=row,ncols=col,figsize=(10,10),dpi=300,sharex=True)
    plt.suptitle("Water transit time distribution - P1 tunnel",fontweight="bold")
    grp_num = 0
    for r in range(row):
        for c in range(col):
            data = []

            current_model = models[c]
            group_name = f"Group {group_number[r]} {current_model}"
            # combine data from all sims
            for sim in range(1,6):
                sim_name = f"{names[r]}A_{current_model}_{sim}"
                _sim_time = combined[sim_name]
                data.extend(_sim_time)

            time_in_ps = [x * 10 for x in data]
            print(f"{group_name} \nrow_col in plot = {r},{c}")
            top_10_slow = sorted(time_in_ps,reverse=True)[:10]
            top_10_fast = sorted(time_in_ps)[:10]

            print(f"TOP 10 frame number of slowest events \n ------------- \n {top_10_slow} \n"
                  f"\n"
                  f"TOP 10 frame number of fastest events \n ------------- \n {top_10_fast} \n")

            slowest=f"Slowest = {max(time_in_ps)}ps"  # Time of the event which took the longest
            fastest = f"Fastest = {min(time_in_ps)}ps"  # Time of the event which is the fastest

            if "opc" in group_name:
                color = overall_color[0]
            elif "tip3p" in group_name:
                color = overall_color[1]
            else:
                color = overall_color[2]

            # https://seaborn.pydata.org/generated/seaborn.histplot.html  <- see here for different 'stat' options
            sns.histplot(data=time_in_ps, ax=ax[r, c], stat="count", binwidth=10, binrange=(0, 400)
                         , color=color, edgecolor='black', log_scale=False,kde=False)

            # axvline is the median of the plot
            # ax[r,c].axvline(sum(time_in_ps)/len(time_in_ps),color='red',linewidth=2,linestyle='dashed')

            ax[r,c].text(0.45,0.9,slowest,transform=ax[r,c].transAxes,va='top')
            # ax[r,c].text(0.65,0.9,fastest,transform=ax[r,c].transAxes,va='center',fontsize=7)

            ax[r,c].set_xlim(0,400)  # x limit capped to 400 for uniform scale
            ax[r,c].set_ylim(0,200)
            ylabel = ax[r,c].get_ylabel()
            ax[r,c].set_ylabel("")
            if ylabel == 'Percent':
                ax [r,c].set_ylim(0,30)
            grp_num += 1

    if ylabel == 'Count':  # statistics used by box histogram to plot the final plot
        y_limits_count_plot = [175,1500,6800,15200]
        for r in range(row):
            for c in range(3):
                lim = (0,y_limits_count_plot[r])
                ax[r,c].set_ylim(lim)

    ax[3, 1].set_xlabel("Number of Events")
    ax[0, 0].set_title("OPC", fontweight='bold')
    ax[0, 1].set_title("TIP3P", fontweight='bold')
    ax[0, 2].set_title("TIP4P-Ew", fontweight='bold')
    ax[0, 0].set_ylabel("TCG1")
    ax[1, 0].set_ylabel("TCG2")
    ax[2, 0].set_ylabel("Time (ps)" + "\n\nTCG3")
    ax[3, 0].set_ylabel("TCG4")
    plt.tight_layout()
    plt.savefig(f"/home/aravind/PhD_local/dean/figures/transit_time/hist_log_trt_{ylabel}manuscript.png")
    print("File saved")

def histogram_matching_frames(tt_results:str,sim_results:str,gd:dict):
    """
    Plot the histogram of frame numbers for events when it's occurring. So, if more than one water
    is entering or releasing from the tunnel, the count increases per histogram bin.
    :param tt_results: Location of Transport tools results
    :param sim_results: Location of simulation results
    :return: saves a plot of histogram in the specified location
    """
    frames = tt_events.get_transit_time(tt_results=tt_results, simulation_results=sim_results,
                                        groups_definitions=gd,
                                        frame_numbers=True)
    combined = frames[0]
    names = list(combined.keys())
    overall_color = sns.color_palette('deep', 3)
    row, col = (4, 3)
    fig,ax = plt.subplots(nrows=row,ncols=col,figsize=(10,8),dpi=300)
    sns.set_context(context="paper",font_scale=1.5)
    plt.suptitle("Matching frames of event occurrences - P1 tunnel",fontweight="bold")
    grp_num = 0
    for r in range(row):
        for c in range(col):
            group_name = names[grp_num]
            data = combined[group_name]
            # time_in_ps = [x * 20 for x in data]
            print(f" {group_name} \nrow_col = {r},{c}")
            # print(time_in_ps)
            if "opc" in group_name:
                color = overall_color[0]
            elif "tip3p" in group_name:
                color = overall_color[1]
            else:
                color = overall_color[2]
            counts = np.bincount(data)
            avg_count = np.mean(counts)
            std_count = np.std(counts)
            # bar plot
            # sns.barplot(x=[''],y=[avg_count],yerr=[std_count],ax=ax[r,c],color=color)
            # hist plot
            plot = sns.histplot(data=data, ax=ax[r, c], stat="count", binwidth=1,binrange=(0, 20000), color=color,
                         edgecolor='black')
            ylabel = ax[r,c].get_ylabel()
            ax[r,c].set_ylabel("")
            grp_num += 1
    ax[3, 1].set_xlabel("Frame Number")
    ax[0, 0].set_title("OPC", fontweight='bold')
    ax[0, 1].set_title("TIP3P", fontweight='bold')
    ax[0, 2].set_title("TIP4P-Ew", fontweight='bold')
    ax[0, 0].set_ylabel("TCG1")
    ax[1, 0].set_ylabel("TCG2")
    ax[2, 0].set_ylabel(ylabel + "\n\n TCG3" )
    ax[3, 0].set_ylabel("TCG4")
    plt.tight_layout()
    plt.savefig(f"/home/aravind/PhD_local/dean/figures/transit_time/hist_frames_{ylabel}_manuscript.png")
    print(frames)


def main():
    tt_results = "/data/aravindramt/dean/tt/tt_0_9_5"
    unassigned_split = "/home/aravind/PhD_local/dean/figures/transport_tools/unassigned_events_sep.csv"
    groups_def = {"P1": [1, 2, 5, 7, 12, 30, 31], "P2": [3, 4, 6, 8, 11, 16, 25, 27, 41, 43, 44, 50, 58],
                  "P3": [10]}
    main_tunnel = {"P1": [1, 2, 5, 7, 12, 30, 31]}
    simulation_results = "/data/aravindramt/dean/md/simulations/"
    sim_results = "/data/aravindramt/dean/md/simulations/"
    save_loc = "/home/aravind/PhD_local/dean/figures/transit_time"

    # Generate data
    # Step1 - Consolidate events - This will generate CSVs of assigned and unassigned events
    tt_events.consolidate_results(tt_results=tt_results, groups_definitions=groups_def,
                           save_location="/home/aravind/PhD_local/dean/figures/transport_tools/")
    #
    # Step2 - Plots
    # PLOT consolidated events results
    plot_consolidated_result('~/PhD_local/dean/figures/transport_tools/consolidated_results.csv',
                             '~/PhD_local/dean/figures/transport_tools/consolidated_unassigned.csv',
                             plot_normalized=False)

    # PLOT stacked entry/release
    # for model in ['opc', 'tip3p', 'tip4pew']:
    #     plot_results_per_tunnel(tt_results=tt_results, groups_definitions=groups_def
    #                             , model=model, unassigned_csv=unassigned_split)

    # Plot transit time
    histogram_count_transit_time(tt_results=tt_results, sim_results=simulation_results, gd=main_tunnel)
    histogram_matching_frames(tt_results=tt_results, sim_results=simulation_results, gd=main_tunnel)

if __name__ == '__main__':
    main()