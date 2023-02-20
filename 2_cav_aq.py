# -*- coding: utf-8 -*-
# Process and plot overall aquaduct and CAVER results

import os
from collections import defaultdict
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def get_caver_data(path_of_input):
    current_dir = path_of_input
    list_of_dirs = os.listdir(path_of_input)
    caver_data = defaultdict(int)
    number_of_tunnels = 0
    for directory in list_of_dirs:
        caver_results_dir = os.path.join(current_dir, (directory + "/caver_analyses/final_clustering/data"))
        try:
            os.chdir(caver_results_dir)
        except NotADirectoryError:
            pass
        except FileNotFoundError:
            pass
        # check for clustering txt file
        if "clustering.txt" in os.listdir(os.getcwd()):
            pass
        else:
            print("No clustering.txt found in" + os.getcwd() + " check the parameters ?")
        with open("clustering.txt") as txt:
            read = txt.readline()
            txt.close()
            epoch = os.getcwd().split("/")
            caver_data[epoch[6]] = int(read.split(' ')[0])
            number_of_tunnels = number_of_tunnels + int(read.split(' ')[0])
    print(f"total tunnels are {number_of_tunnels}")
    caver_sorted = [(k, v) for k, v in caver_data.items()]
    # Sort by 1A, 1.4A, 1.8A, 2.4A and 3A .. OPC, TIP3P, TIP4P-Ew
    caver_sorted.sort(key=lambda x: (x[0].split("_")[1], x))
    opc = caver_sorted[0:25]
    opc.sort(key=lambda x: (x[0].split("_")[0][:-1], x))
    tip3p = caver_sorted[25:50]
    tip3p.sort(key=lambda x: (x[0].split("_")[0][:-1], x))
    tip4p = caver_sorted[50:75]
    tip4p.sort(key=lambda x: (x[0].split("_")[0][:-1], x))
    opc_out = dict(opc)
    tip3p_out = dict(tip3p)
    tip4p_out = dict(tip4p)

    def split_five(input):
        result = []
        temp_list = []
        counter = 1
        for key, value in input.items():
            temp_list.append(value)
            if counter % 5 == 0:
                result.append(temp_list)
                temp_list = []
            counter += 1
        if temp_list:
            result.append(temp_list)
        return result

    opc = split_five(opc_out)
    tip3p = split_five(tip3p_out)
    tip4p = split_five(tip4p_out)
    return opc, tip3p, tip4p


def read_aquaduct_input():
    """
    Read 5_analysis_results.txt get 'Number of traceable residues' for the MD results (OLD FILE PATTERN)
    """
    models = ['opc', 'tip3p', 'tip4pew']
    epoch_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
    sim_list = ["1", "2", "3", "4", "5"]
    water_count = 0
    model_wat = []
    for model in models:
        epoch_wat = []
        model_water_count=0
        for epoch in epoch_list:
            epoch_water_count = 0
            sim_wat = []
            for sim in sim_list:
                aq_results_dir = '/mnt/NAS/dean_water_models/md/' + model + '/' + epoch + '/' + sim + '/aquaduct/'
                with open(aq_results_dir + "5_analysis_results.txt") as txt:
                    read = txt.readlines()
                    # print(read)
                    txt.close()
                    line_number = read.index("Names of traced molecules: WAT\n")
                    sim_wat.append(int(read[line_number + 2].split(sep=':')[1]))
                    water_count += int(read[line_number + 2].split(sep=':')[1])
                    epoch_water_count += int(read[line_number + 2].split(sep=':')[1])
                    model_water_count += int(read[line_number + 2].split(sep=':')[1])
            epoch_wat.append([*sim_wat])
            print(f"{epoch} - {epoch_water_count}")
        model_wat.append([*epoch_wat])
        print(f"{model}- {model_water_count}")
    print(f"Overall waters are {water_count}")
    return model_wat


def plot_aquaduct(water_data):
    """
    Consolidated plot of groups + models + waters traced
    :param water_data: The output from read_aquaduct_input()
    :return: Plots
    """
    group_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
    width = 0.6
    fig, ax = plt.subplots(1, 5, figsize=(12, 3), dpi=300)
    plt.subplots_adjust(top=0.8, left=0.07, right=0.89, wspace=0.5, bottom=0.2)
    fig.suptitle('AQUA-DUCT WATERS TRACED', fontsize=12, fontweight='bold')
    x = np.arange(1, 6)
    for epoch in range(5):
        ax[epoch].set_title(r'{}Å'.format(group_list[epoch][:-1]), size=12)
        ax[epoch].set_xticks(x)
        ax[epoch].yaxis.grid(linestyle='--', color='gray', alpha=0.7, zorder=-1)
        opc_x = x - width / 3
        opc_height = water_data[0][epoch]
        tip3p_x = x
        tip3p_height = water_data[1][epoch]
        tip4pew_x = x + width / 3
        tip4pew_height = water_data[2][epoch]
        ax[epoch].bar(opc_x, opc_height, width / 3, label='OPC')
        ax[epoch].bar(tip3p_x, tip3p_height, width / 3, label='TIP3P')
        ax[epoch].bar(tip4pew_x, tip4pew_height, width / 3, label='TIP4P-Ew')
        ax[epoch].tick_params(axis='both', which='major', labelsize=10)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10)
    ax[2].set_xlabel('SIMULATION ID', fontsize=10)
    ax[0].set_ylabel('NUMBER OF WATERS', fontsize=10)
    plt.savefig("/home/aravind/PhD_local/dean/figures/caver&aquaduct/aquaduct.png")
    # plt.show()


def plot_aquaduct_per_group(water_data):
    """
    Subplots of waters per group per simulation
    :param water_data: The output from read_aquaduct_input()
    :return:
    """
    epoch_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
    plt.style.use("grayscale")
    fig, ax = plt.subplots(5, 3, figsize=(8.27, 11.7), dpi=300)
    opc = water_data[0]
    tip3p = water_data[1]
    tip4pew = water_data[2]

    def make_df(listoflist):
        dic = {"1A": [], "1.4A": [], "1.8A": [], "2.4A": [], "3A": []}
        keys = list(dic.keys())
        i = 0
        for sim in listoflist:
            dic[keys[i]] = sim
            i += 1
        return pd.DataFrame(dic, index=[1, 2, 3, 4, 5])

    opc_df = make_df(opc)
    tip3p_df = make_df(tip3p)
    tip4pew_df = make_df(tip4pew)
    # plt.subplots_adjust(wspace=0.4,hspace=0.3)
    plt.suptitle('AQUA-DUCT Traced Waters', fontsize=15, fontweight='bold')
    # fig.text(s='Sim Number',x = 0.48, y=0.07,fontsize=15,fontweight='bold')
    # fig.text(s='Number of waters', x=0.05, y=0.37, fontsize=15,rotation=90,fontweight='bold')
    opc_df.plot(ax=ax[0, 0], kind='bar', y='1A', legend=False, rot=0,
                color='r').set_ylabel("Group 1A", fontsize=15)
    opc_df.plot(ax=ax[1, 0], kind='bar', y='1.4A', legend=False, rot=0,
                color='r').set_ylabel("Group 1.4A", fontsize=15)
    opc_df.plot(ax=ax[2, 0], kind='bar', y='1.8A', legend=False,
                color='r').set_ylabel("Group 1.8A", fontsize=15)
    opc_df.plot(ax=ax[3, 0], kind='bar', y='2.4A', legend=False, rot=0,
                color='r').set_ylabel("Group 2.4A", fontsize=15)
    opc_df.plot(ax=ax[4, 0], kind='bar', y='3A', legend=False, rot=0,
                color='r').set_ylabel("Group 3A", fontsize=15)
    tip3p_df.plot(ax=ax[0, 1], kind='bar', y='1A', legend=False, rot=0, color='g')
    tip3p_df.plot(ax=ax[1, 1], kind='bar', y='1.4A', legend=False, rot=0, color='g')
    tip3p_df.plot(ax=ax[2, 1], kind='bar', y='1.8A', legend=False, rot=0, color='g')
    tip3p_df.plot(ax=ax[3, 1], kind='bar', y='2.4A', legend=False, rot=0, color='g')
    tip3p_df.plot(ax=ax[4, 1], kind='bar', y='3A', legend=False, rot=0, color='g')
    tip4pew_df.plot(ax=ax[0, 2], kind='bar', y='1A', legend=False, rot=0, color='c')
    tip4pew_df.plot(ax=ax[1, 2], kind='bar', y='1.4A', legend=False, rot=0, color='c')
    tip4pew_df.plot(ax=ax[2, 2], kind='bar', y='1.8A', legend=False, rot=0, color='c')
    tip4pew_df.plot(ax=ax[3, 2], kind='bar', y='2.4A', legend=False, rot=0, color='c')
    tip4pew_df.plot(ax=ax[4, 2], kind='bar', y='3A', legend=False, rot=0, color='c')
    ax[0, 0].set_title("P1", fontsize=15, fontweight='bold')
    ax[0, 1].set_title("P2", fontsize=15, fontweight='bold')
    ax[0, 2].set_title("P3", fontsize=15, fontweight='bold')
    ax[4, 1].set_xlabel("Simulation Number", fontsize=10, fontweight='bold')
    for axes in ax.flatten():
        axes.tick_params(which='both', size=3, length=5)
    plt.tight_layout(pad=1, w_pad=0.5, h_pad=0.9)
    plt.savefig("/home/aravind/PhD_local/dean/figures/caver&aquaduct/waters_per_group.png")


def plot_caver(OPC_caver, TIP3P_caver, TIP4Pew_caver):
    epoch_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
    width = 0.6
    fig, ax = plt.subplots(1, 5, figsize=(12, 3), dpi=300)
    plt.subplots_adjust(top=0.8, left=0.07, right=0.89, wspace=0.5, bottom=0.2)
    fig.suptitle('CAVER TUNNELS - OVERALL PER GROUP', fontsize=12, fontweight='bold')
    x = np.arange(1, 6)
    for epoch in range(5):
        ax[epoch].set_title(r'{}Å'.format(epoch_list[epoch][:-1]), size=12)
        ax[epoch].set_xticks(x)
        ax[epoch].yaxis.grid(linestyle='--', color='gray', alpha=0.7, zorder=-1)
        opc_x = x - width / 3
        opc_height = OPC_caver[epoch]
        tip3p_x = x
        tip3p_height = TIP3P_caver[epoch]
        tip4pew_x = x + width / 3
        tip4pew_height = TIP4Pew_caver[epoch]
        ax[epoch].bar(opc_x, opc_height, width / 3, label='OPC')
        ax[epoch].bar(tip3p_x, tip3p_height, width / 3, label='TIP3P')
        ax[epoch].bar(tip4pew_x, tip4pew_height, width / 3, label='TIP4P-Ew')
        ax[epoch].tick_params(axis='both', which='major', labelsize=10)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10)
    ax[2].set_xlabel('SIMULATION ID', fontsize=10)
    ax[0].set_ylabel('TOTAL CAVER TUNNELS', fontsize=10)
    plt.savefig("/home/aravind/PhD_local/dean/figures/caver&aquaduct/caver.png")
    # plt.show()


if __name__ == '__main__':
    sim_dir = "/data/aravindramt/dean/md/simulations"

    # AQUA-DUCT
    # ------------------------------
    waters = read_aquaduct_input()
    # plot_aquaduct_per_group(waters)
    plot_aquaduct(waters)

    # CAVER
    # ------------------------------
    opc, tip3p, tip4pew = get_caver_data(sim_dir)
    plot_caver(opc, tip3p, tip4pew)
