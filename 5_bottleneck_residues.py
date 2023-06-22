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
    bottlenecks = defaultdict(dict)
    # bottlenecks = []
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
                bottlenecks[sc_id] = filtered_freq
                # bottlenecks=bottlenecks+list(filtered_freq.keys())
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

def process_comparative_resutls(comparative_results:dict):
    all_residues = set()
    opc= defaultdict(int)
    tip3p = defaultdict(int)
    tip4pew = defaultdict(int)
    for group in comparative_results:
        group_values = comparative_results[group]
        for id in group_values:
            id_values: dict = group_values[id]
            for res,freq in id_values.items():
                all_residues.add(res)
                number_of_frames = int(20000 * freq)
                if 'opc' in group:
                    if res in opc:
                        opc[res] += number_of_frames
                    else:
                        opc[res] = number_of_frames
                elif 'tip3p' in group:
                    if res in tip3p:
                        tip3p[res] += number_of_frames
                    else:
                        tip3p[res] = number_of_frames
                elif 'tip4p' in group:
                    if res in tip4pew:
                        tip4pew[res] += number_of_frames
                    else:
                        tip4pew[res] = number_of_frames

    # organize
    for residue in all_residues:
        if residue not in opc:
            opc[residue] = 0
        if residue not in tip3p:
            tip3p[residue] = 0
        if residue not in tip4pew:
            tip4pew[residue] = 0

    def convert2df(object: dict) -> pd.DataFrame:
        df = pd.DataFrame.from_dict(object, orient='index', columns=['Frames'])
        df.index = pd.to_numeric(df.index)
        df = df.sort_index()
        filtered_df = df[df['Frames'] > 1000]
        return filtered_df

    opc_df = convert2df(opc)
    tip3p_df = convert2df(tip3p)
    tip4pew_df = convert2df(tip4pew)
    return [opc_df,tip3p_df,tip4pew_df]

def plot_comparative(plot_data,save_location:str) -> None:
    sns.set_context(context='paper',font_scale=0.8)
    colors = sns.color_palette('deep', 3)
    opc_df,tip3p_df,tip4pew_df = plot_data
    row,col =1,3
    fig, axes = plt.subplots(nrows=row, ncols=col, dpi=300)
    sns.barplot(ax=axes[0],data=opc_df,y=opc_df.index,x='Frames',orient='h',color=colors[0])
    sns.barplot(ax=axes[1],data=tip3p_df, y=tip3p_df.index, x='Frames', orient='h',color=colors[1])
    sns.barplot(ax=axes[2],data=tip4pew_df, y=tip4pew_df.index, x='Frames', orient='h',color=colors[2])
    axes[0].set_ylabel("Residue Number")
    axes[0].set_title("OPC")
    axes[1].set_title("TIP3P")
    axes[2].set_title("TIP4P-Ew")
    plt.suptitle("Bottleneck residues - P1 tunnel")
    plt.tight_layout()
    plt.savefig(save_location+"p1.png")


if __name__ == '__main__':

    result_file = "/data/aravindramt/dean/tt/tt_0_9_5_bottleneck/statistics/" \
                  "2-filtered_tunnels_statistics_bottleneck_residues.txt"
    p1 = [1, 2, 5, 7, 12, 30, 31]
    p2 = [3, 4, 6, 8, 11, 16, 25, 27, 41, 43, 44, 50, 58]
    p3 = [10]

    # frames by model
    save_location_comparative = "/home/aravind/PhD_local/dean/figures/bottlenecks/"
    comparative_analysis_results = "/data/aravindramt/dean/tt/tt_0_9_5_bottleneck/statistics/comparative_analysis"
    result = get_comparative_results(comparative_analysis_results,p1)
    plot_data= process_comparative_resutls(result)
    plot_comparative(plot_data,save_location_comparative)


