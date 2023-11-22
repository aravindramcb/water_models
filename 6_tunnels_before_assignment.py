# -*- coding: utf-8 -*-
# Python script to plot the presence of tunnels for given groups of SuperCluster IDs
import os
import pickle
from collections import defaultdict
import pandas
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter


def _process_csv(csv_file):
    """
    Reads super_cluster CSV files and returns sim_id and number of frames
    Generally the CSV file only has the simIDs which are participating in the formation of that supercluster but, for
    plotting the overall figure, I also need to have the other simulations which are not a member of that superclsuter
    so is function will return all the sim_ids which are participating and not participating for easier plotting.
    :param csv_file: CSV profile file
    :return:dict of md_tag:number_of_frames
    """
    simulations = ('1A_tip3p_1', '3A_opc_1', '1.4A_tip3p_2', '2.4A_tip3p_4', '3A_tip4pew_2', '2.4A_tip3p_3',
                   '1.4A_tip3p_5', '1A_tip4pew_3', '2.4A_tip3p_5', '1A_tip3p_5', '1.4A_tip3p_1', '3A_tip4pew_3',
                   '2.4A_opc_2', '2.4A_opc_1', '3A_opc_2', '2.4A_tip4pew_3', '1.4A_opc_2', '3A_tip3p_2', '1A_opc_5',
                   '3A_tip3p_3', '1A_tip4pew_4', '1A_opc_4', '1.4A_tip3p_3', '1.8A_tip3p_1', '3A_tip3p_4', '3A_tip3p_5',
                   '1.4A_tip4pew_3', '2.4A_tip3p_1', '2.4A_opc_3', '3A_tip4pew_4', '2.4A_tip4pew_4', '3A_tip3p_1',
                   '1A_tip4pew_5', '1A_tip4pew_1', '1A_tip3p_2', '2.4A_opc_4', '1.4A_tip4pew_5', '1.4A_opc_1',
                   '1.4A_tip3p_4', '1.4A_tip4pew_1', '1.8A_opc_5', '1.8A_tip4pew_1', '2.4A_tip3p_2', '1.8A_tip3p_5',
                   '2.4A_tip4pew_5', '1.8A_tip3p_3', '1A_opc_1', '1.4A_tip4pew_4', '3A_tip4pew_1', '1A_tip3p_3',
                   '1A_tip3p_4', '1.4A_opc_5', '1.8A_tip4pew_5', '1A_opc_3', '3A_tip4pew_5', '2.4A_tip4pew_1',
                   '3A_opc_3', '1.4A_opc_4', '1.8A_opc_4', '1A_opc_2', '3A_opc_4', '1.8A_tip4pew_3', '2.4A_tip4pew_2',
                   '1.8A_opc_1', '1.8A_tip3p_2', '2.4A_opc_5', '1.4A_tip4pew_2', '1.8A_tip4pew_4', '3A_opc_5',
                   '1.8A_opc_2', '1A_tip4pew_2', '1.8A_tip3p_4', '1.4A_opc_3', '1.8A_tip4pew_2', '1.8A_opc_3')

    df = pd.read_csv(csv_file, header=0, index_col=None, usecols=[0, 1])
    df.reset_index(inplace=True)
    df = df.iloc[:, [0, 1]].rename(columns={df.columns[0]: 'Md_traj', df.columns[1]: 'snapshot'})
    # Get MD_label and unique snapshots
    md_label_frames = defaultdict(set)
    for index, values in df.iterrows():
        md_label = values.Md_traj
        snapshot = values.snapshot
        md_label_frames[md_label].add(snapshot)
    # Count snapshots and return count
    number_frames = defaultdict(int)
    for sim_id in simulations:
        if sim_id in md_label_frames:
            frame_count = len(md_label_frames[sim_id])
            number_frames[sim_id] = frame_count
        else:
            number_frames[sim_id] = 0
    return number_frames


def _get_others_scids(tunnels_dict: dict, tt_results_location: str):
    """
    Find SCIDs of the IDs which are not in tunnels_dict
    :param tunnels_dict: Tunnels dictionary
    :param tt_results_location: TransportTools results location
    :param save_location: Save location of dump file
    :return: Dump file of the other tunnels
    """
    csv_file_location = os.path.join(tt_results_location, "data", "super_clusters", "CSV_profiles", "filtered01")
    file_names = os.listdir(csv_file_location)
    sc_ids = {int(filename.split('_')[2].split('.')[0]) for filename in file_names}
    sc_ids_in_tunnels_dict = {id for sublist in tunnels_dict.values() for id in sublist}
    other_scids = sc_ids - sc_ids_in_tunnels_dict
    return other_scids


def count_frames(tunnels_dict, tt_results_location, save: bool = False, save_location: str = None):
    """
    Counts number of frames per simulation ID for the given supercluster ID and consolidates it per tunnel (P1,P2,P3)
    :param tunnels_dict: Tunnel definitions
    :return: dataframe of tunnels and their frame numbers
    """
    from tqdm import tqdm
    csv_file_location = os.path.join(tt_results_location, "data", "super_clusters", "CSV_profiles", "filtered01")
    other_ids = _get_others_scids(tunnels_dict, tt_results_location)
    tunnels_dict['others'] = other_ids
    tunnels = defaultdict(pandas.DataFrame)
    for tunnel in tunnels_dict:
        scids = tunnels_dict[tunnel]
        df_scid = pd.DataFrame()
        for scid in tqdm(scids, desc=f"Processing tunnel {tunnel}"):
            # for scid in scids:
            if scid < 10:
                tunnel_id = f"0{scid}"
                csv_filename = f"super_cluster_{tunnel_id}.csv"
                csv_file = os.path.join(csv_file_location, csv_filename)
            else:
                csv_filename = f"super_cluster_{scid}.csv"
                csv_file = os.path.join(csv_file_location, csv_filename)
            num_frames = _process_csv(csv_file)
            num_frames_df = pd.DataFrame.from_dict(num_frames, orient='index', columns=[scid])
            df_scid = pd.concat([df_scid, num_frames_df], axis=1)
        ids = list(df_scid.index)
        ids.sort(key=lambda x: (x.split("_")[0][:-1], x))
        df_scid = df_scid.reindex(ids)
        tunnels[tunnel] = df_scid
    if save:
        # Save file because it is a time-consuming operation
        save_file = os.path.join(save_location, 'number_frames.pkl')
        with open(save_file, 'wb') as f:
            pickle.dump(tunnels, f)
    return tunnels


def plot(data: dict = None, dump_file: str = None, save_location: str = None):
    """
    Plot the tunnels per group + others
    :param data: dictionary containing the number of frames
    :param save_location: save locaton
    :return: none
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    if dump_file:
        with open(dump_file, 'rb') as f:
            data: dict = pickle.load(f)
    elif data:
        pass
    else:
        raise FileNotFoundError(f"{dump_file} not found")
    colors = sns.color_palette('deep', 3)
    row, col = (5, 4)
    fig, ax = plt.subplots(nrows=row, ncols=col, figsize=(11.7, 8.9), dpi=300, sharex='col')
    sns.set(style='whitegrid')
    sns.set_palette('deep')
    plt.suptitle("TUNNEL OCCURRENCES BEFORE ASSIGNMENT", fontsize=12, fontweight='bold')
    tunnels = ["P1", "P2", "P3", "others"]
    for _col in range(col):
        i, j, k = 0, 5, 10
        for _row in range(row):
            df = data[tunnels[_col]]
            df_sum = df.sum(axis=1)
            opc = df_sum[i:i + 5]
            opc.reset_index(drop=True, inplace=True)
            tip3p = df_sum[j:j + 5]
            tip3p.reset_index(drop=True, inplace=True)
            tip4pew = df_sum[k:k + 5]
            tip4pew.reset_index(drop=True, inplace=True)
            plot_df = pd.concat([opc, tip3p, tip4pew], axis=1)
            print(plot_df)
            sns.barplot(data=plot_df, ax=ax[_row, _col], width=0.5, errorbar='se', capsize=.1, linewidth=1, errwidth=1)
            if _col == 0:  # P1
                ax[_row,_col].set_ylim(0, 25000)
            elif _col == 1:  # P2
                ax[_row, _col].set_ylim(0, 25000)
            elif _col == 2:  # P3
                ax[_row, _col].set_ylim(0, 1250)
            elif _col == 3:  # Others
                ax[_row, _col].set_ylim(0, 12500)
            ax[_row,_col].set_xticklabels("")
            i += 15
            j += 15
            k += 15

    ax[4,2].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    # ax[4,2].set_ylim(0,15)
    ax[0, 0].set_title("P1", fontsize=10, fontweight='bold')
    ax[0, 1].set_title("P2", fontsize=10, fontweight='bold')
    ax[0, 2].set_title("P3", fontsize=10, fontweight='bold')
    ax[0, 0].set_ylabel("TCG0", fontsize=10)
    ax[1, 0].set_ylabel("TCG1", fontsize=10)
    ax[2, 0].set_ylabel("TCG2", fontsize=10)
    ax[3, 0].set_ylabel("TCG3", fontsize=10)
    ax[4, 0].set_ylabel("TCG4", fontsize=10)
    ax[0, 3].set_title("Others", fontsize=10, fontweight='bold')
    plt.tight_layout(pad=1.9, w_pad=0.3, h_pad=0.3)
    opc_patch = Patch(color=colors[0], label='OPC')
    tip3p_patch = Patch(color=colors[1], label='TIP3P')
    tip4pew_patch = Patch(color=colors[2], label='TIP4P-Ew')
    ax[0, 3].legend([opc_patch, tip3p_patch, tip4pew_patch], ["OPC", "TIP3P", "TIP4P-Ew"],
                    bbox_to_anchor=(-0.55, 1.05, 0, 0), loc="lower right",
                    borderaxespad=1, ncol=3)
    fig.text(x=0.01, y=0.4, s="NUMBER OF SNAPSHOTS\n\n", rotation=90)
    fig.text(x=0.4, y=0.01, s="TUNNELS & MODELS")
    plt.savefig(save_location + "before.png")


def main():
    tt_results_location = "/data/aravindramt/dean/tt/tt_0_9_5"
    save_location = '/home/aravind/PhD_local/dean/figures/transport_tools/'
    groups_def = {"P1": {1, 2, 5, 7, 12, 30, 31}, "P2": {3, 4, 6, 8, 11, 16, 25, 27, 41, 43, 44, 50, 58},
                  "P3": {10}}
    # count_frames can be used with or without saving dump file, here I have saved a dump because my CSV files are
    # very large (>4GB), so I cannot wait every time to modify results.
    # count_frames(tunnels_dict=groups_def, tt_results_location=tt_results_location, save=True,
    #              save_location=save_location)
    dump_file_name=os.path.join(save_location,'number_frames.pkl')
    plot(dump_file=dump_file_name,save_location=save_location)


if __name__ == '__main__':
    main()
