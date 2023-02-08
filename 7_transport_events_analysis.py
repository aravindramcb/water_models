# Python script to process the results of TransportTools and generate plots for Overall results and comparative analysis
# -*- coding: utf-8 -*-
__author__ = 'Aravind Selvaram Thirunavukarasu'
__email__ = 'arathi@amu.edu.pl, aravind1233@gmail.com'

import itertools
from collections import defaultdict
import seaborn as sns
import pandas as pd
import os
from dataclasses import dataclass, field
from matplotlib import pyplot as plt


@dataclass(order=True, repr=True, unsafe_hash=False)
class TTEventsStats:
    SC_ID: int = field(repr=True, hash=True)
    No_Sims: int = field(repr=True, hash=True)
    Total_No_Frames: int = field(repr=True, hash=True)
    Avg_No_Frames: float = field(repr=True, hash=True)
    Avg_BR: float = field(repr=True, hash=True)
    StDev_BR: float = field(repr=True, hash=True)
    Max_BR: float = field(repr=True, hash=True)
    Avg_Len: float = field(repr=True, hash=True)
    StDev_Len: float = field(repr=True, hash=True)
    Avg_Cur: float = field(repr=True, hash=True)
    StDev_Cur: float = field(repr=True, hash=True)
    Avg_throug: float = field(repr=True, hash=True)
    StDev_through: float = field(repr=True, hash=True)
    Priority: float = field(repr=True, hash=True)
    Num_Events: int = field(repr=True, hash=True)
    Num_entries: int = field(repr=True, hash=True)
    Num_releases: int = field(repr=True, hash=True)


@dataclass(order=True, repr=True, unsafe_hash=False)
class tt_unassigned:
    Sim_ID: str
    Num_Events: int
    Num_entries: int
    Num_releases: int


@dataclass(order=True, repr=True, unsafe_hash=False)
class before_assignment:
    SC_ID: int
    No_Sims: int
    Total_No_Frames: int
    Avg_No_Frames: float
    Avg_BR: float
    StDev_BR: float
    Max_BR: float
    Avg_Len: float
    StDev_Len: float
    Avg_Cur: float
    StDev_Cur: float
    Avg_Throug: float
    StDev_Throug: float
    Priority: float



def _find_sc_per_group(comparative_results_dir: str):
    """
    Reads the comparative analysis results and then groups according to user definitions. Also gives out
    'others' - which are the SCs not in any groups and 'unassigned' which are not assigned by transport tools

    :param comparative_results_dir: Directory in which comparative analysis results are present
    :return:dict(comparative_group_name:[events_assigned_to_scs,unassigned_events]
    """
    # Parse the results and store it
    comparative_analysis_results = [d for d in os.listdir(comparative_results_dir) if
                                    os.path.isdir(os.path.join(comparative_results_dir, d))]
    comparative_analysis_results.sort()
    parsed_results = defaultdict(list)
    for group in comparative_analysis_results:
        file_path = os.path.join(comparative_results_dir, group, "4-filtered_events_statistics.txt")

        if os.path.isfile(file_path):
            with open(file_path, "r") as file:
                lines = file.readlines()
                events_per_group = []
                for line in lines[21:-1]:  # start reading from line 22
                    # Parse the line and create a DataClass object
                    if "-" in line:
                        continue
                    else:
                        data = [float(x) if '.' in x else int(x) for x in line.split(",")]
                        events_per_group.append(TTEventsStats(*data))
                        # print the obj
                # print(events_per_group[0])
                unasigned = lines[-1].split(",")
                total = int(unasigned[0].split()[5])
                entry = int(unasigned[1])
                release = int(unasigned[2])
                group_name = group
                unasigned_per_group = tt_unassigned(group_name, total, entry, release)
                # print(unasigned_per_group)
                parsed_results[group_name] = [events_per_group, unasigned_per_group]
        else:
            print(f"{file_path} does not exist.")
    return parsed_results


def _process_unassigned(comparative_results_dir: str):
    """
    Reads the comparative analysis results and then groups according to user definitions.
    :param comparative_results_dir: Directory in which comparative analysis results are present
    :return:dict
    """
    # Parse the results and store it
    comparative_analysis_results = [d for d in os.listdir(comparative_results_dir) if
                                    os.path.isdir(os.path.join(comparative_results_dir, d))]
    comparative_analysis_results.sort()
    parsed_results = defaultdict(list)
    for group in comparative_analysis_results:
        file_path = os.path.join(comparative_results_dir, group, "2-filtered_tunnels_statistics.txt")

        if os.path.isfile(file_path):
            with open(file_path, "r") as file:
                lines = file.readlines()
                events_per_group = []
                for line in lines[18:-1]:  # start reading from line 22
                    # Parse the line and create a DataClass object
                    if "-" in line:
                        continue
                    else:
                        data = [float(x) if '.' in x else int(x) for x in line.split(",")]
                        events_per_group.append(before_assignment(*data))
                group_name = group
                parsed_results[group_name] = events_per_group
        else:
            print(f"{file_path} does not exist.")
    return parsed_results


def _get_scids_of_groups(comparative_analysis_results, groups_definition: dict = None, show_info: bool = False):
    """
    Gets the SC_IDs of user defined groups if they are present in comparative analysis results.
    :param show_info: Print the SC_IDs selected per user defined group.
    :param comparative_analysis_results: results dir location
    :param groups_definition: User defined definition
    :return: SC_IDs of tunnels split by user defined groups + others + unassigned
    """
    # Split by user defined groups definitions
    parsed_results = _find_sc_per_group(comparative_analysis_results)
    sc_id_by_tunnels = defaultdict(list)
    group_names = list(groups_definition.keys())
    for epoch in parsed_results.keys():
        assigned = parsed_results[epoch][0]  # 0 -> assigned , 1 -> unassigned
        SC_IDs = [x.SC_ID for x in assigned]
        # Split into groups, which are defined at the beginning of script
        for user_group_name in group_names:  # for P1 in [p1,p2,p3]
            sc_id_in_group = [x for x in SC_IDs if x in groups_definition[user_group_name]]
            if epoch not in sc_id_by_tunnels:
                # Create a new entry in sc_id_by_tunnels if comparative_group does not exist
                sc_id_by_tunnels[epoch] = [{user_group_name: sc_id_in_group}]
            else:
                # If it exists, append other values in user group . eg., P2, P3.. others
                sc_id_by_tunnels[epoch][0][user_group_name] = sc_id_in_group
        # Add a new entry called 'others' where it has all the SCIDs which does not belong to the groups which user
        # defined
        selected_sc_ids = [item for sublist in sc_id_by_tunnels[epoch][0].values() for item in sublist]
        other_scids = [x for x in SC_IDs if x not in selected_sc_ids]
        sc_id_by_tunnels[epoch][0]["others"] = other_scids

        # Add 'unassigned' to the results
        unassigned = parsed_results[epoch][1]

        # reorder
        order = group_names + ["others"]
        sc_id_by_tunnels[epoch][0] = {key: sc_id_by_tunnels[epoch][0][key] for key in order
                                      if key in sc_id_by_tunnels[epoch][0]}
    if show_info:
        for group in sc_id_by_tunnels:
            # print(f"for group {group}")
            print(sc_id_by_tunnels[group])
            # print(" ")

    return sc_id_by_tunnels


def _get_unassigned_events(tt_results: str):
    """
    Returns only unassigned events per comparative groups definition / folder names
    :param tt_results: TransportTools results directory.
    :return:dict{comparative_group_name):unassigned_values,.....}
    """
    outlier_file = os.path.join(tt_results, "data", "super_clusters", "details", "outlier_transport_events_details.txt")
    with open(outlier_file, 'r') as infile:
        contents = infile.readlines()
    unassigned = defaultdict(dict)
    current_category = 'entry'
    for line in contents[4:]:
        # entry: (from Simulation: AQUA-DUCT ID, (Resname:Residue), start_frame->end_frame; ... )
        # 'from 3A_opc_3: 1952, (WAT:8178), 19098->19131;'
        if 'from' in line and not 'release' in line:
            line = line.strip()
            sim_ID = line.split(" ")[1].split(":")[0]
            num_events = len(line.split(";")[:-1])
            if current_category in unassigned[sim_ID]:
                unassigned[sim_ID][current_category].append(num_events)
            else:
                unassigned[sim_ID] = {"entry": [0], "release": [0]}
                unassigned[sim_ID][current_category].append(num_events)
        elif 'release' in line:
            current_category = "release"
    return unassigned


def _get_data_from_TT(tt_super_cluster_details):
    """
    Helper function to parse the details from filtered_supercluster_details1.txt file and returns the entry and release
    len() for every sim #.
    :param tt_super_cluster_details:
    The file from tt_results/data/super_clusters
    /details/filtered_super_clusters_details1.txt
    :returns super_cluster_IDS , dict(entry_water), dict(release_water)
    """
    super_cluster_ids = []
    entry_water = defaultdict(list)
    release_water = defaultdict(list)
    caver_ids = []
    with open(tt_super_cluster_details, 'r') as results_file:
        _read = results_file.readlines()
        readfile = [i.strip("\n") for i in _read]
        parsed_data = {}
        i = 0
        while i < len(readfile) - 1:
            _super_cluster_id = []  # Temp super cluster id
            _epoch_water_entry = defaultdict(int)  # Temp water entry
            _epoch_water_release = defaultdict(int)  # Temp water release
            while not readfile[i].startswith('Super'):
                i += 1
            if readfile[i].startswith('Super'):
                sc_id = int(readfile[i].split(" ")[2])
                _super_cluster_id.append(sc_id)  # Add to tmp list to pass below
                super_cluster_ids.append(sc_id)  # Add to master list

                # print(readfile[i])
            while not readfile[i].startswith('entry'):  # Move the pointer to entry:
                if readfile[i].startswith('release'):  # if there is no entry, move to release
                    break
                i += 1
            while not readfile[i].startswith('release'):  # Move pointer until release, and parse the entry data

                if readfile[i].startswith('entry'):  # overriding for mismatch of entries in the file
                    pass
                else:
                    if readfile[i].startswith("-"):  # break loop if there is no entry
                        break
                    epoch, waters = readfile[i].split(":", 1)
                    sim_id = epoch.split(sep=" ")[1]  # omit "from "
                    entry_waters = waters.split(sep=";")
                    # print(entry_waters)
                    entry_water_count = len(entry_waters) - 1
                    _epoch_water_entry[sim_id] = entry_water_count
                    # print("ENTRY-"+sim_id+"="+str(entry_water_count))
                i += 1

            while not readfile[i].startswith('-'):  # Move the pointer until the dash part and parse the release data
                if readfile[i].startswith('Super'):
                    break
                if readfile[i].startswith('release'):
                    pass
                    # print(readfile[i])
                else:
                    # print("release", i)
                    epoch, waters = readfile[i].split(":", 1)
                    sim_id = epoch.split(sep=" ")[1]
                    release_waters = waters.split(sep=";")
                    # print(release_waters)
                    release_water_count = len(release_waters) - 1
                    _epoch_water_release[sim_id] = release_water_count
                    # print("RELEASE-"+sim_id+"="+str(release_water_count))
                #     One cycle of entry and release gets over here, i+1 to goto next supercluster ID
                i += 1
            entry_water[sc_id].append(_epoch_water_entry)
            release_water[sc_id].append(_epoch_water_release)
    results_file.close()
    return super_cluster_ids, entry_water, release_water


def _process_results_for_SCs(required_SCIDs: list, tt_results: str, model: str):
    """
    Helper function
    Give me the required Supercluster IDs and I will give you the input and release number of events for them.
    :param tt_results: Location of TransportTools results
    :param model: Water model "OPC" or "TIP3P" or "TIP4PEW" (in small letters)
    :param required_SCIDs: Super Cluster IDs to be processed. Eg., [1] or [1,3,5] etc.,
    """
    tt_sc_details2_file = os.path.join(tt_results, "data", "super_clusters", "details",
                                       "filtered_super_cluster_details2.txt")
    sim_list = ['1A_opc_1', '1A_opc_2', '1A_opc_3', '1A_opc_4', '1A_opc_5', '1.4A_opc_1', '1.4A_opc_2', '1.4A_opc_3',
                '1.4A_opc_4', '1.4A_opc_5', '1.8A_opc_1', '1.8A_opc_2', '1.8A_opc_3', '1.8A_opc_4',
                '1.8A_opc_5', '2.4A_opc_1', '2.4A_opc_2', '2.4A_opc_3', '2.4A_opc_4', '2.4A_opc_5', '3A_opc_1',
                '3A_opc_2', '3A_opc_3', '3A_opc_4', '3A_opc_5',

                '1A_tip3p_1', '1A_tip3p_2', '1A_tip3p_3', '1A_tip3p_4', '1A_tip3p_5', '1.4A_tip3p_1',
                '1.4A_tip3p_2', '1.4A_tip3p_3', '1.4A_tip3p_4', '1.4A_tip3p_5', '1.8A_tip3p_1', '1.8A_tip3p_2',
                '1.8A_tip3p_3', '1.8A_tip3p_4', '1.8A_tip3p_5', '2.4A_tip3p_1', '2.4A_tip3p_2',
                '2.4A_tip3p_3', '2.4A_tip3p_4', '2.4A_tip3p_5', '3A_tip3p_1', '3A_tip3p_2', '3A_tip3p_3', '3A_tip3p_4',
                '3A_tip3p_5',

                '1A_tip4pew_1', '1A_tip4pew_2', '1A_tip4pew_3', '1A_tip4pew_4', '1A_tip4pew_5',
                '1.4A_tip4pew_1', '1.4A_tip4pew_2', '1.4A_tip4pew_3', '1.4A_tip4pew_4', '1.4A_tip4pew_5',
                '1.8A_tip4pew_1', '1.8A_tip4pew_2', '1.8A_tip4pew_3', '1.8A_tip4pew_4', '1.8A_tip4pew_5',
                '2.4A_tip4pew_1', '2.4A_tip4pew_2', '2.4A_tip4pew_3', '2.4A_tip4pew_4', '2.4A_tip4pew_5',
                '3A_tip4pew_1', '3A_tip4pew_2', '3A_tip4pew_3', '3A_tip4pew_4', '3A_tip4pew_5']

    if model == "opc":
        sim_list = sim_list[0:25]
    elif model == "tip3p":
        sim_list = sim_list[25:50]
    elif model == "tip4pew":
        sim_list = sim_list[50:75]

    # print(model, sim_list)
    sc_ids, entry, release = _get_data_from_TT(tt_sc_details2_file)
    _tmp_df = pd.DataFrame(0, index=[0], columns=sim_list)
    # ENTRY
    for sc_ID in required_SCIDs:
        # location = list(entry.keys()).index(sc_ID)
        #  Add 'entry' column to the dataframe
        sc_id_df = pd.DataFrame(entry[sc_ID])
        _tmp_df = pd.concat([_tmp_df.convert_dtypes(), sc_id_df.convert_dtypes()], axis=0, ignore_index=True)
    # _tmp_df --> Dataframe of Transport events per SC, can be used to plot the highest contributor
    # print(_tmp_df)
    sum_entry_per_replica = _tmp_df.sum(axis=0)
    tt_entry_df_for_all_sc = pd.DataFrame(sum_entry_per_replica, columns=['Entry_events'])
    _tmp_df = pd.DataFrame(0, index=[0], columns=sim_list)
    # RELEASE
    for sc_ID in required_SCIDs:
        # location = list(entry.keys()).index(sc_ID)
        sc_id_df = pd.DataFrame(release[sc_ID])
        _tmp_df = pd.concat([_tmp_df.convert_dtypes(), sc_id_df.convert_dtypes()], axis=0, ignore_index=True)
        # print(_tmp_df)
    # _tmp_df --> Dataframe of Transport events per SC, can be used to plot the highest contributor
    sum_release = _tmp_df.sum(axis=0)
    tt_release_df_for_all_sc = pd.DataFrame(sum_release, columns=['Release_events'])
    return tt_entry_df_for_all_sc, tt_release_df_for_all_sc


def plot_tunnels_before_assignment(groups_definitons: dict, comparative_analysis_results_folder,
                                   save_location: str = None):
    # Get the Total_No_Frames and normalize it
    parsed_results = _process_unassigned(comparative_analysis_results_folder)
    frames_per_group = defaultdict(dict)
    scids_of_gd = list(groups_definitons.values())
    full_values = list(itertools.chain(*scids_of_gd))
    for group_name in parsed_results:
        frames = defaultdict(dict)
        current_group = parsed_results[group_name]  # Current group = opc_1
        others = []
        for gd in groups_definitons:  # For P1 in P1,P2,P3
            frames_in_group = []

            # Iterate through each scid in current group and add it to list if present, else to others.
            for event in current_group:
                scid = event.SC_ID
                total_frames = event.Total_No_Frames
                if scid in groups_definitons[gd]:
                    frames_in_group.append(total_frames)
                elif scid not in full_values:
                    others.append(total_frames)
            sum_per_gd = int(sum(frames_in_group))
            sum_others_per_gd = int(sum(others))
            # frames[gd]=sum_per_gd
            frames.update({gd: sum_per_gd})
        frames.update(dict(others=int(sum(others))))
        frames_per_group[group_name] = frames
    df = pd.DataFrame.from_dict(frames_per_group)

    # Plot output
    # ------------
    import matplotlib.patches as mpatches
    colors = sns.color_palette('deep', 3)
    row, col = (5, 4)
    fig, ax = plt.subplots(nrows=5, ncols=4, figsize=(8.27, 11.7), dpi=300, sharex='col')
    sns.set(style='whitegrid')
    sns.set_palette('deep')
    plt.suptitle("TUNNEL OCCURRENCES BEFORE ASSIGNMENT", fontsize=12, fontweight='bold')
    o, t, p = (0, 5, 10)
    for r in range(row):
        for c in range(col):
            plot_data = df.iloc[c:c + 1, [o, t, p]]
            sns.barplot(data=plot_data, ax=ax[r, c])
            ax[r, c].set_xticklabels([])
            # set y limits
            # if c ==0:
            #     ax[r,c].set_ylim(0,24000)
            # elif c == 1:
            #     ax[r,c].set_ylim(0,17000)  #17000
            # elif c == 2:
            #     ax[r,c].set_ylim(0,6000)  #6000
            # elif c ==3:
            #     ax[r,c].set_ylim(0,30000)  #30000
        o += 1
        t += 1
        p += 1

    ax[0, 0].set_title("P1", fontsize=10, fontweight='bold')
    ax[0, 1].set_title("P2", fontsize=10, fontweight='bold')
    ax[0, 2].set_title("P3", fontsize=10, fontweight='bold')
    ax[0, 0].set_ylabel("Group 1A", fontsize=10)
    ax[1, 0].set_ylabel("Group 1.4A", fontsize=10)
    ax[2, 0].set_ylabel("Group 1.8A", fontsize=10)
    ax[3, 0].set_ylabel("Group 2.4A", fontsize=10)
    ax[4, 0].set_ylabel("Group 3A", fontsize=10)
    ax[0, 3].set_title("Others", fontsize=10, fontweight='bold')
    plt.tight_layout(pad=1.9, w_pad=0.3, h_pad=0.3)
    opc_patch = mpatches.Patch(color=colors[0], label='OPC')
    tip3p_patch = mpatches.Patch(color=colors[1], label='TIP3P')
    tip4pew_patch = mpatches.Patch(color=colors[2], label='TIP4P-Ew')
    ax[0, 3].legend([opc_patch, tip3p_patch, tip4pew_patch], ["OPC", "TIP3P", "TIP4P-Ew"],
                    bbox_to_anchor=(-0.55, 1.05, 0, 0), loc="lower right",
                    borderaxespad=1, ncol=3)
    fig.text(x=0.01, y=0.5, s="NUMBER OF SNAPSHOTS", rotation=90)
    fig.text(x=0.4, y=0.01, s="TUNNELS & MODELS")
    save_file = os.path.join(save_location, "tunnel_occurrence_no_scale.png")
    plt.savefig(save_file)


def consolidate_results(tt_results: str, groups_definitions: dict, save_location=None):
    """
    Process the TT assigned and unassigned events and process results to produce entry,release and consolidated results.
    :param tt_results: TransportTools results location
    :param groups_definitions: The definitions of groups formed by SCs.
    :returns consolidated_df Dataframe
    [optional] saving of consolidated df to consolidated_results.csv
    """
    global sc_ids
    sim_list = ['1A_opc_1', '1A_opc_2', '1A_opc_3', '1A_opc_4', '1A_opc_5', '1.4A_opc_1', '1.4A_opc_2', '1.4A_opc_3',
                '1.4A_opc_4', '1.4A_opc_5', '1.8A_opc_1', '1.8A_opc_2', '1.8A_opc_3', '1.8A_opc_4',
                '1.8A_opc_5', '2.4A_opc_1', '2.4A_opc_2', '2.4A_opc_3', '2.4A_opc_4', '2.4A_opc_5', '3A_opc_1',
                '3A_opc_2', '3A_opc_3', '3A_opc_4', '3A_opc_5',

                '1A_tip3p_1', '1A_tip3p_2', '1A_tip3p_3', '1A_tip3p_4', '1A_tip3p_5', '1.4A_tip3p_1',
                '1.4A_tip3p_2', '1.4A_tip3p_3', '1.4A_tip3p_4', '1.4A_tip3p_5', '1.8A_tip3p_1', '1.8A_tip3p_2',
                '1.8A_tip3p_3', '1.8A_tip3p_4', '1.8A_tip3p_5', '2.4A_tip3p_1', '2.4A_tip3p_2',
                '2.4A_tip3p_3', '2.4A_tip3p_4', '2.4A_tip3p_5', '3A_tip3p_1', '3A_tip3p_2', '3A_tip3p_3', '3A_tip3p_4',
                '3A_tip3p_5',

                '1A_tip4pew_1', '1A_tip4pew_2', '1A_tip4pew_3', '1A_tip4pew_4', '1A_tip4pew_5',
                '1.4A_tip4pew_1', '1.4A_tip4pew_2', '1.4A_tip4pew_3', '1.4A_tip4pew_4', '1.4A_tip4pew_5',
                '1.8A_tip4pew_1', '1.8A_tip4pew_2', '1.8A_tip4pew_3', '1.8A_tip4pew_4', '1.8A_tip4pew_5',
                '2.4A_tip4pew_1', '2.4A_tip4pew_2', '2.4A_tip4pew_3', '2.4A_tip4pew_4', '2.4A_tip4pew_5',
                '3A_tip4pew_1', '3A_tip4pew_2', '3A_tip4pew_3', '3A_tip4pew_4', '3A_tip4pew_5']

    comparative_analysis_results = os.path.join(tt_results, "statistics", "comparative_analysis")
    consolidated_df = pd.DataFrame()
    rows, cols = (5, 4)
    # tunnel_id = ["P1", "P2", "P3"]
    tunnel_id = list(groups_definitions.keys()) + ["others"]
    group = [1, 2, 3, 4, 5]

    # CONSOLIDATE ENTRY + RELEASE PER SIMULATION
    groups_scs = _get_scids_of_groups(comparative_analysis_results,
                                      groups_definitions)  # get corresponding SC per group
    opc_contents = list(groups_scs.keys())[0:5]  # Assuming folder 1-5 is OPC and so on
    tip3p_contents = list(groups_scs.keys())[5:10]
    tip4pew_contents = list(groups_scs.keys())[10:15]
    for model in ['opc', 'tip3p', 'tip4pew']:
        i, j = (0, 5)
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
                # print("TUNNEL ID =", tunnel_id[col])
                entry_df, release_df = _process_results_for_SCs(sc_ids, tt_results, model=model)
                entry_data = entry_df[i:j]  # To plot entry per tunnel
                release_data = release_df[i:j]
                _combined_df = pd.concat([entry_data, release_data], axis=1)
                # _combined_df.reset_index(inplace=True)
                print(_combined_df.sum())
                sum_list = list(_combined_df.sum(axis=1))
                print(sum_list)
                col_name = tunnel_id[col] + "_" + model + "_" + str(group[(row)])
                _tmp_df = pd.DataFrame({col_name: sum_list})
                consolidated_df = pd.concat([consolidated_df, _tmp_df], axis=1)
            i += 5
            j += 5

    # Process Unassigned TT events

    unassigned_df = pd.DataFrame()
    unassigned = _get_unassigned_events(tt_results=tt_results)
    # '1.4A_opc_1': {'entry': 1, 'release': 0}
    for sim in sim_list:
        if sim in unassigned:
            num_entry = unassigned[sim]['entry']
            num_release = unassigned[sim]['release']
            _tmp_unas_df = pd.DataFrame({sim: [sum(num_entry), sum(num_release)]}, index=['Entry', 'Release'])
            unassigned_df = pd.concat([unassigned_df, _tmp_unas_df], axis=1)
        else:
            _tmp_unas_df = pd.DataFrame({sim: [0, 0]}, index=['Entry', 'Release'])
            unassigned_df = pd.concat([unassigned_df, _tmp_unas_df], axis=1)
    unassigned_df.to_csv(save_location + 'unassigned_events_sep.csv')

    # Consolidate unassigned and group by group and model
    _consolidated_unassigned = pd.DataFrame()
    sum_unassigned_df = unassigned_df.sum(axis=0)
    unassigned_grouped_by_group = pd.DataFrame()
    models = ['opc','tip3p','tip4pew']
    current_model = 0
    group_number = 0
    for k in range(0, len(sum_unassigned_df), 5):
        chunk = sum_unassigned_df[k:k + 5]
        chunk = chunk.reset_index(drop=True)
        chunk = chunk.to_frame()
        chunk.columns = [f'{models[current_model]}_{group_number+1}']
        unassigned_grouped_by_group = pd.concat([unassigned_grouped_by_group, chunk], axis=1)
        if k == 20:
            current_model +=1
        elif k ==45:
            current_model +=1
        group_number += 1
        if group_number == 5:
            group_number -= 5
    unassigned_grouped_by_group.to_csv(save_location+'consolidated_unassigned.csv',index=False)
    if save_location is not None:
        consolidated_df.to_csv(save_location + 'consolidated_results.csv', index=False)

    return consolidated_df


def plot_consolidated_result(consolidated_csv_file: str, unassigned_csv: str, plot_normalized: bool = False):
    import matplotlib.patches as mpatches
    colors = sns.color_palette('deep', 3)
    consolidated_df = pd.read_csv(consolidated_csv_file)
    i, j, k = (0, 20, 40)  #
    rows, cols = (5, 4)
    import matplotlib.pyplot as plt
    sns.set(style='white')
    fig, ax = plt.subplots(nrows=5, ncols=5, figsize=(12.27, 11.7), dpi=150, sharex='col')
    if plot_normalized is True:
        plt.suptitle("Events Scalled by TIP3P & Consolidated by tunnels and models", fontsize=15, fontweight='bold')
    else:
        plt.suptitle("Events Consolidated by tunnels and models", fontsize=15, fontweight='bold')
    for row in range(rows):
        for col in range(cols):
            _df = consolidated_df.iloc[:, [i, j, k]]
            if plot_normalized == True:
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
        sns.barplot(data=_unas_df, ax=ax[unas_row, 4], width=0.5)
        ax[unas_row, 4].set_xticklabels([])
        a += 1
        b += 1
        c += 1
    plt.tight_layout(pad=3, w_pad=0.1, h_pad=0.1)
    opc_patch = mpatches.Patch(color=colors[0], label='OPC')
    tip3p_patch = mpatches.Patch(color=colors[1], label='TIP3P')
    tip4pew_patch = mpatches.Patch(color=colors[2], label='TIP4P-Ew')
    ax[0, 4].legend([opc_patch, tip3p_patch, tip4pew_patch], ["OPC", "TIP3P", "TIP4P-Ew"],
                    bbox_to_anchor=(-0.9, 1.15, 0, 0), loc="lower right",
                    borderaxespad=1, ncol=3)
    ax[0, 0].set_title("P1", fontsize=10, fontweight='bold')
    ax[0, 1].set_title("P2", fontsize=10, fontweight='bold')
    ax[0, 2].set_title("P3", fontsize=10, fontweight='bold')
    ax[0, 0].set_ylabel("Group 1A", fontsize=10)
    ax[1, 0].set_ylabel("Group 1.4A", fontsize=10)
    ax[2, 0].set_ylabel("Group 1.8A", fontsize=10)
    ax[3, 0].set_ylabel("Group 2.4A", fontsize=10)
    ax[4, 0].set_ylabel("Group 3A", fontsize=10)
    ax[0, 3].set_title("Others", fontsize=10, fontweight='bold')
    ax[0, 4].set_title("Unassigned", fontsize=10, fontweight='bold')
    fig.text(x=0.01, y=0.5, s="NUMBER OF SNAPSHOTS", rotation=90)
    fig.text(x=0.45, y=0.01, s="TUNNELS & MODELS")

    if plot_normalized == True:
        plt.savefig('/home/aravind/PhD_local/dean/figures/transport_tools/consolidated_tt_events_scaled.png')
    else:
        plt.savefig('/home/aravind/PhD_local/dean/figures/transport_tools/consolidated_tt_mean_se_events.png')
    # plt.show()


def plot_results_per_tunnel(comparative_analysis_results: str, groups_definitions: dict, model: str,
                            tt_sc_details2_file_loc: str, unassigned_csv: str, plot_type: str = "consolidated"):
    global sc_ids, index
    import matplotlib.pyplot as plt
    sns.set(style='white')
    colors = {"opc": ['xkcd:dusk', 'xkcd:cool blue'], "tip3p": ['xkcd:dirt', 'xkcd:browny orange'],
              "tip4pew": ['xkcd:moss green', 'xkcd:seaweed']}
    fig, ax = plt.subplots(nrows=5, ncols=5, figsize=(10.27, 11.7), dpi=150, sharex='col')
    plt.suptitle(f"{model}".upper() + " Events", fontsize=15, fontweight='bold')
    rows, cols = (5, 4)
    i, j = (0, 5)
    x_label = [1, 2, 3, 4, 5]
    groups_scs = _get_scids_of_groups(comparative_analysis_results,
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
            entry_df, release_df = _process_results_for_SCs(sc_ids, tt_sc_details2_file_loc, model)
            data = entry_df[i:j]  # To plot entry per tunnel
            release_data = release_df[i:j]
            _combined_df = pd.concat([data, release_data], axis=1)
            _combined_df.reset_index(inplace=True, drop=True)
            _combined_df["Total"] = _combined_df.sum(axis=1)

            # using seaborn to plot
            sns.barplot(ax=ax[row, col], data=_combined_df, x=_combined_df.index, y='Total', color=colors[model][0])
            sns.barplot(ax=ax[row, col], data=_combined_df, x=_combined_df.index, y='Release_events',
                        color=colors[model][1])

            # CONSOLIDATED PLOT per model
            # plot_consolidated = sns.barplot(data=_combined_df.sum(axis=1), ax=ax[row, col])
            # TO plot entry and release events, un-comment and modify below
            # plot = sns.barplot(data=data,x=data.index,y='Entry_events',ax=ax[row,col],color='b')

            # # Set x and y ticks
            ax[row, col].set_ylabel("")
            ax[row, col].set_xlabel("")
            ax[row, col].set_xticks([0, 1, 2, 3, 4])
            ax[row, col].set_xticklabels(x_label, rotation=0)

        i += 5
        j += 5

    if model == 'opc':
        index = 1
    elif model == 'tip3p':
        index = 6
    elif model == 'tip4pew':
        index = 11

    # Load unassigned data
    unassigned_df = pd.read_csv(unassigned_csv)
    unassigned_df.rename(index={0: 'entry', 1: 'exit'})

    # Plot unassigned
    for unas_col in range(5):
        _unas_df = unassigned_df.iloc[:, [index]]
        _unas_df = _unas_df.T
        _unas_df = _unas_df.rename(columns={0: "entry", 1: "release"})
        _unas_df['total'] = _unas_df.sum(axis=1)
        _unas_df.reset_index(inplace=True)
        sns.barplot(data=_unas_df, x='index', y='total', ax=ax[unas_col, 4], width=0.2, color=colors[model][0])
        sns.barplot(data=_unas_df, x='index', y='release', ax=ax[unas_col, 4], width=0.2, color=colors[model][1])
        ax[unas_col, 4].set_ylabel("")
        ax[unas_col, 4].set_xlabel("")
        ax[unas_col, 4].set_xticks([0])
        ax[unas_col, 4].set_xticklabels(["Unassigned"])
        index += 1

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
    ax[0, 0].set_ylabel("Group 1A", fontsize=15)
    ax[1, 0].set_ylabel("Group 1.4A", fontsize=15)
    ax[2, 0].set_ylabel("Group 1.8A", fontsize=15)
    ax[3, 0].set_ylabel("Group 2.4A", fontsize=15)
    ax[4, 0].set_ylabel("Group 3A", fontsize=15)
    ax[4, 2].set_xlabel("Simulation #", fontsize=15)

    # if model == "tip3p":
    #     ax[1,2].text(0.5, 0.5,'No events', ha='center', va='center', transform=ax[1,2].transAxes)
    #     ax[1, 3].text(0.5, 0.5, 'No events', ha='center', va='center', transform=ax[1, 3].transAxes)

    plt.savefig(f'/home/aravind/PhD_local/dean/figures/transport_tools/combined_events_{model}.png')
    # plt.savefig(f'/home/aravind/PhD_local/dean/figures/{model}_consolidated.png')
    # plt.show()


def main():
    tt_results = "/mnt/gpu/dean/tt/tt_0_9_5"
    unasigned = "/home/aravind/PhD_local/dean/figures/transport_tools/unassigned_individual.csv"

    # groups_def = {"P1": [1, 2, 5, 7, 12, 30, 31], "P2": [3, 4, 6, 11, 25, 27, 41, 43, 44, 50, 58, 16],
    #               "P3": [8, 9, 10, 24]}
    groups_def = {"P1": [1, 2, 5, 7, 12, 30, 31], "P2": [3, 4, 6, 8, 9, 11, 16, 24, 25, 27, 41, 43, 44, 50, 58],
                  "P3": [10
                      , 29, 37, 47, 69, 78]}
    # DEBUG
    # _get_data_from_TT(results)
    # _process_results_for_SCs([1],results,'opc')
    # _get_scids_of_groups(comparative_results, groups_def,show_info=True)

    # Tunnels before transport events assignment - from comparative analysis results
    # plot_tunnels_before_assignment(groups_def,comparative_results,'/home/aravind/PhD_local/dean/figures/transport_tools/')

    # Transport events analysis
    # Step1 - Consolidate results - This will generate CSVs of assigned and unassigned events
    consolidate_results(tt_results=tt_results, groups_definitions=groups_def,
                        save_location="/home/aravind/PhD_local/dean/figures/transport_tools/")

    # Step2 - Plots
    # PLOT consolidated results
    # plot_consolidated_result('~/PhD_local/dean/figures/transport_tools/consolidated_results_with_unassigned.csv',
    #                          '~/PhD_local/dean/figures/transport_tools/unassigned_events.csv',
    #                          plot_normalized=False)

    # PLOT stacked entry/release
    # plot_results_per_tunnel(comparative_analysis_results=comparative_results, groups_definitions=groups_def
    #                         , model="tip4pew", tt_sc_details2_file_loc=results, unassigned_csv=unasigned)


if __name__ == '__main__':
    main()
