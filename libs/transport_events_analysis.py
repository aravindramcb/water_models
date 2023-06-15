# Python script to process the results of TransportTools
# -*- coding: utf-8 -*-
__author__ = 'Aravind Selvaram Thirunavukarasu'
__email__ = 'arathi@amu.edu.pl, aravind1233@gmail.com'

from collections import defaultdict
import pandas as pd
import os
from dataclasses import dataclass, field
from statistics import mean


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
class TTUnassigned:
    Sim_ID: str
    Num_Events: int
    Num_entries: int
    Num_releases: int


@dataclass(order=True, repr=True, unsafe_hash=False)
class BeforeAssignment:
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
                unasigned_per_group = TTUnassigned(group_name, total, entry, release)
                # print(unasigned_per_group)
                parsed_results[group_name] = [events_per_group, unasigned_per_group]
        else:
            print(f"{file_path} does not exist.")
    return parsed_results


def get_scids_of_groups(comparative_analysis_results, groups_definition: dict = None, show_info: bool = False):
    """
    Gets the SC_IDs of user defined groups only if they are present in comparative analysis results.
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
        sc_id = [x.SC_ID for x in assigned]
        # Split into groups, which are defined at the beginning of script
        for user_group_name in group_names:  # for P1 in [p1,p2,p3]
            sc_id_in_group = [x for x in sc_id if x in groups_definition[user_group_name]]
            if epoch not in sc_id_by_tunnels:
                # Create a new entry in sc_id_by_tunnels if comparative_group does not exist
                sc_id_by_tunnels[epoch] = [{user_group_name: sc_id_in_group}]
            else:
                # If it exists, append other values in user group . eg., P2, P3.. others
                sc_id_by_tunnels[epoch][0][user_group_name] = sc_id_in_group
        # Add a new entry called 'others' where it has all the SCIDs which does not belong to the groups which user
        # defined
        selected_sc_ids = [item for sublist in sc_id_by_tunnels[epoch][0].values() for item in sublist]
        other_scids = [x for x in sc_id if x not in selected_sc_ids]
        sc_id_by_tunnels[epoch][0]["others"] = other_scids
        # reorder
        order = group_names + ["others"]
        sc_id_by_tunnels[epoch][0] = {key: sc_id_by_tunnels[epoch][0][key] for key in order
                                      if key in sc_id_by_tunnels[epoch][0]}
    if show_info:
        for group in sc_id_by_tunnels:
            print(f"for group {group}")
            print(sc_id_by_tunnels[group])
            # print(" ")

    return sc_id_by_tunnels


def _get_unassigned_events(tt_results: str):
    """
    Returns only unassigned events per comparative groups definition / folder names
    :param tt_results: TransportTools results directory.
    :return:dict{comparative_group_name:unassigned_values,.....}
    """
    outlier_file = os.path.join(tt_results, "data", "super_clusters", "details", "outlier_transport_events_details.txt")
    with open(outlier_file, 'r') as infile:
        contents = infile.readlines()
    unassigned = defaultdict(dict)
    current_category = 'entry'
    for line in contents[4:]:
        # entry: (from Simulation: AQUA-DUCT ID, (Resname:Residue), start_frame->end_frame; ... )
        # 'from 3A_opc_3: 1952, (WAT:8178), 19098->19131;'
        if 'from' in line:
            if 'release' not in line:
                line = line.strip()
                sim_id = line.split(" ")[1].split(":")[0]
                num_events = len(line.split(";")[:-1])
                if current_category in unassigned[sim_id]:
                    unassigned[sim_id][current_category] += num_events
                else:
                    unassigned[sim_id] = {"entry": 0, "release": 0}
                    unassigned[sim_id][current_category] += num_events
            elif 'release' in line:
                current_category = "release"
        elif 'release' in line:
            current_category = "release"
    return unassigned


def _get_data_from_tt(tt_super_cluster_details, debug: bool = False, frame_numbers: bool = False):
    """
    Helper function to parse the details from filtered_supercluster_details2.txt file and returns the entry and release
    len() for every sim #.
    :param tt_super_cluster_details:
    The file from tt_results/data/super_clusters
    /details/filtered_super_clusters_details1.txt
    :returns dict(entry,release) of Aquaduct frames per event , dict(entry_water), dict(release_water)
    """
    # EVENTS
    num_frames_aq_water = defaultdict(list)
    frame_numbers_aq_water = defaultdict(list)
    entry_water = defaultdict(list)
    release_water = defaultdict(list)
    if debug:
        print("---------------------------------------------------------------")
    with open(tt_super_cluster_details, 'r') as results_file:
        _read = results_file.readlines()
        readfile = [i.strip("\n") for i in _read]
        i = 0
        while i < len(readfile) - 1:
            _super_cluster_id = []  # Temp super cluster id
            _epoch_water_entry = defaultdict(int)  # Temp water entry
            _epoch_water_release = defaultdict(int)  # Temp water release
            _epoch_frames_entry = defaultdict(list)
            _epoch_frames_release = defaultdict(list)
            _epoch_frame_numbers_entry = defaultdict(list)
            _epoch_frame_numbers_release = defaultdict(list)

            while not readfile[i].startswith('Super'):
                i += 1

            if readfile[i].startswith('Super'):
                sc_id = int(readfile[i].split(" ")[2])
                _super_cluster_id.append(sc_id)  # Add to tmp list to pass below
                # print(readfile[i])

            while not readfile[i].startswith('from'):  # move pointer below "Tunnel clusters:"
                i += 1

            _tunnel_clusters_sim_list = []
            while not readfile[i] == '':
                sim_id = readfile[i].split()[1][:-1]
                _tunnel_clusters_sim_list.append(sim_id)
                i += 1

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
                    # print(f"Parsing data for SCID - {sc_id}")
                    frames = [i.split(sep=',')[2] for i in entry_waters[:-1]]
                    number_of_frames_entry = []
                    frame_numbers_entry = []
                    for element in frames:
                        start, end = element.split("->")
                        number_of_frames = int(end) - int(start)
                        if frame_numbers:
                            _frame_numbers = list(range(int(start), int(end) + 1))
                            frame_numbers_entry = frame_numbers_entry + _frame_numbers  # Frame 1,2,3,4 etc.
                        number_of_frames_entry.append(number_of_frames)  # count of frames
                    entry_water_count = len(entry_waters) - 1
                    if sim_id in _tunnel_clusters_sim_list:
                        _epoch_water_entry[sim_id] = entry_water_count
                        _epoch_frames_entry[sim_id] = number_of_frames_entry
                        if frame_numbers:
                            _epoch_frame_numbers_entry[sim_id] = frame_numbers_entry

                    if debug:
                        # print("There is no tunnel cluster for the entry event below :  ")
                        print(f"[ENTRY] - SCID {sc_id} - {sim_id} - {entry_water_count}")
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
                    frames = [i.split(sep=',')[2] for i in release_waters[:-1]]
                    number_of_frames_release = []
                    frame_numbers_release = []
                    for element in frames:
                        start, end = element.split("->")
                        number_of_frames = int(end) - int(start)
                        if frame_numbers:
                            _frame_numbers = list(range(int(start), int(end) + 1))
                            frame_numbers_release = frame_numbers_release + _frame_numbers
                        number_of_frames_release.append(number_of_frames)
                    release_water_count = len(release_waters) - 1
                    if sim_id in _tunnel_clusters_sim_list:
                        _epoch_frames_release[sim_id] = number_of_frames_release
                        _epoch_water_release[sim_id] = release_water_count
                        if frame_numbers:
                            _epoch_frame_numbers_release[sim_id] = frame_numbers_release
                    if debug:
                        print("There is no tunnel cluster for the release event below :  ")
                        print(f"[RELEASE] - SCID {sc_id} - {sim_id} - {release_water_count}")
                    # print("RELEASE-"+sim_id+"="+str(release_water_count))
                #     One cycle of entry and release gets over here, i+1 to goto next supercluster ID
                i += 1
            entry_water[sc_id].append(_epoch_water_entry)
            release_water[sc_id].append(_epoch_water_release)
            num_frames_aq_water[sc_id] = [_epoch_frames_entry, _epoch_frames_release]
            if frame_numbers:
                frame_numbers_aq_water[sc_id] = [_epoch_frame_numbers_entry, _epoch_frame_numbers_release]
    results_file.close()
    if debug:
        print("----------------------------------------------------------------\n \n")
    if frame_numbers:
        return frame_numbers_aq_water
    else:
        return num_frames_aq_water, entry_water, release_water


def get_entry_release_events(required_SCIDs: list, tt_results: str, model: str):
    """
    Helper function
    Give me the required Supercluster IDs and I will give you the entry and release number of events for them.
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
    super_cluster_ids, entry, release = _get_data_from_tt(tt_sc_details2_file)
    _tmp_df = pd.DataFrame(0, index=[0], columns=sim_list)
    # ENTRY
    for sc_ID in required_SCIDs:
        # location = list(entry.keys()).index(sc_ID)
        #  Add 'entry' column to the dataframe
        sc_id_df = pd.DataFrame(entry[sc_ID])
        try:
            _tmp_df = pd.concat([_tmp_df.convert_dtypes(), sc_id_df.convert_dtypes()], axis=0, ignore_index=True)
        except ValueError:
            pass
    # _tmp_df --> Dataframe of Transport events per SC, can be used to plot the highest contributor
    # print(_tmp_df)
    sum_entry_per_replica = _tmp_df.sum(axis=0)
    tt_entry_df_for_all_sc = pd.DataFrame(sum_entry_per_replica, columns=['Entry_events'])
    _tmp_df = pd.DataFrame(0, index=[0], columns=sim_list)
    # RELEASE
    for sc_ID in required_SCIDs:
        # location = list(entry.keys()).index(sc_ID)
        try:
            sc_id_df = pd.DataFrame(release[sc_ID])
            _tmp_df = pd.concat([_tmp_df.convert_dtypes(), sc_id_df.convert_dtypes()], axis=0, ignore_index=True)
        except ValueError:
            pass
        # print(_tmp_df)
    # _tmp_df --> Dataframe of Transport events per SC, can be used to plot the highest contributor
    sum_release = _tmp_df.sum(axis=0)
    tt_release_df_for_all_sc = pd.DataFrame(sum_release, columns=['Release_events'])
    return tt_entry_df_for_all_sc, tt_release_df_for_all_sc


def consolidate_results(tt_results: str, groups_definitions: dict, save_location=None):
    """
    Process the TT assigned and unassigned events and process results to produce entry,release and consolidated results.
    :param save_location: Location to save
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
    groups_scs = get_scids_of_groups(comparative_analysis_results,
                                     groups_definitions)  # get corresponding SC per group which are water transporters

    opc_folders = list(groups_scs.keys())[0:5]  # Assuming folder 1-5 is OPC and so on
    tip3p_folders = list(groups_scs.keys())[5:10]
    tip4pew_folders = list(groups_scs.keys())[10:15]
    for model in ['opc', 'tip3p', 'tip4pew']:
        i, j = (0, 5)  # parse data from 5 simulations
        for row in range(rows):
            for col in range(cols):
                if model == "opc":
                    group_name = opc_folders[row]
                    print(f"For group {group_name} {tunnel_id[col]}")
                    sc_ids = groups_scs[group_name][0][tunnel_id[col]]
                if model == "tip3p":
                    group_name = tip3p_folders[row]
                    print(f"For group {group_name}")
                    sc_ids = groups_scs[group_name][0][tunnel_id[col]]
                if model == "tip4pew":
                    group_name = tip4pew_folders[row]
                    print(f"For group {group_name}")
                    sc_ids = groups_scs[group_name][0][tunnel_id[col]]
                # print("TUNNEL ID =", tunnel_id[col])
                entry_df, release_df = get_entry_release_events(sc_ids, tt_results, model=model)
                entry_data = entry_df[i:j]  # To plot entry per tunnel
                release_data = release_df[i:j]
                try:
                    _combined_df = pd.concat([entry_data, release_data], axis=1)
                except ValueError:
                    pass
                # _combined_df.reset_index(inplace=True)
                print(_combined_df.sum())
                sum_list = list(_combined_df.sum(axis=1))
                # print(sum_list)
                col_name = tunnel_id[col] + "_" + model + "_" + str(group[row])
                _tmp_df = pd.DataFrame({col_name: sum_list})
                consolidated_df = pd.concat([consolidated_df, _tmp_df], axis=1)
            i += 5
            j += 5

    # Process Unassigned TT events

    unassigned_df = pd.DataFrame()
    unassigned = _get_unassigned_events(tt_results=tt_results)
    # '1.4A_opc_1': {'entry': 1, 'release': 0}
    # Fix missing values - use entry,release as 0,0 if there are no unassigned events in that simulation
    for sim in sim_list:
        if sim in unassigned:
            num_entry = unassigned[sim]['entry']
            num_release = unassigned[sim]['release']
            _tmp_unas_df = pd.DataFrame({sim: [num_entry, num_release]}, index=['Entry', 'Release'])
            unassigned_df = pd.concat([unassigned_df, _tmp_unas_df], axis=1)
        else:
            _tmp_unas_df = pd.DataFrame({sim: [0, 0]}, index=['Entry', 'Release'])
            unassigned_df = pd.concat([unassigned_df, _tmp_unas_df], axis=1)
    unassigned_df.to_csv(save_location + 'unassigned_events_sep.csv')

    # Consolidate unassigned and group by group and model
    _consolidated_unassigned = pd.DataFrame()
    sum_unassigned_df = unassigned_df.sum(axis=0)
    unassigned_grouped_by_group = pd.DataFrame()
    models = ['opc', 'tip3p', 'tip4pew']
    current_model = 0
    group_number = 0
    # iterate through 75 folders
    for k in range(0, len(sum_unassigned_df), 5):
        chunk = sum_unassigned_df[k:k + 5]
        chunk = chunk.reset_index(drop=True)
        chunk = chunk.to_frame()
        chunk.columns = [f'{models[current_model]}_{group_number + 1}']
        unassigned_grouped_by_group = pd.concat([unassigned_grouped_by_group, chunk], axis=1)
        # index 0-20 is OPC 20-45 is TIP3P, 45-75 is TIP4P-Ew
        if k == 20:
            current_model += 1
        elif k == 45:
            current_model += 1
        group_number += 1
        if group_number == 5:
            group_number -= 5
    unassigned_grouped_by_group.to_csv(save_location + 'consolidated_unassigned.csv', index=False)
    if save_location is not None:
        consolidated_df.to_csv(save_location + 'consolidated_results.csv', index=False)

    return consolidated_df


def get_transit_time(tt_results, simulation_results: str, groups_definitions: dict, frame_numbers: bool = False,
                     debug: bool = False):
    """
    Transit time = Time spent inside the tunnel for a single aquaduct event.
    This is to compare how much frames an event takes time to pass through a tunnel.
    :param tt_results Folder in which tt results are present
    :param simulation_results folder in which simulation results are present
    :param groups_definitions Tunnels definition, tunnel_name:superclusters
    :param frame_numbers Return the frame numbers per simulaiton
    :returns Three dictionaries of combined_events, entry, release which contains list(list) of number of frames the
    event took place/ if frame_numbers=True, then returns the range of frame numbers which events took place

    """
    # Read filtered_super_cluster_details_2.txt
    sc_details_file = os.path.join(tt_results, "data", "super_clusters", "details",
                                   "filtered_super_cluster_details2.txt")
    if frame_numbers:
        frames = _get_data_from_tt(sc_details_file, frame_numbers=True)
    else:
        parsed_results = _get_data_from_tt(sc_details_file)
        frames = parsed_results[0]
    directories = [d for d in os.listdir(simulation_results) if os.path.isdir(os.path.join(simulation_results, d))]
    directories.sort(key=lambda x: (x.split("_")[0][:-1], x))
    combined_dict = defaultdict(list)
    entry_dict = defaultdict(list)
    release_dict = defaultdict(list)
    comparative_results_loc = os.path.join(tt_results, "statistics/comparative_analysis")
    name_of_tunnel = list(groups_definitions.keys())[0]
    scids = groups_definitions[name_of_tunnel]
    scids_in_group = get_scids_of_groups(comparative_analysis_results=comparative_results_loc,
                                         groups_definition=groups_definitions, show_info=True)
    for scid in scids:
        current_values = frames[scid]  # current_values[0]=entry [1]=release
        # combine entry and release to a single value
        combined_values = defaultdict(list)
        for d in current_values:
            for key, value in d.items():
                if not len(value) == 0:
                    combined_values[key].extend(value)
        # include simulations with empty events
        combined_values_all = {k: combined_values.get(k, 0) for k in directories}
        index = 0

        # names are the suffixes of the group names, my full group names are -['opc_1', 'opc_1.4', 'opc_1.8', 'opc_2.4',
        # 'opc_3', 'tip3p_1', 'tip3p_1.4', 'tip3p_1.8', 'tip3p_2.4', 'tip3p_3', 'tip4pew_1', 'tip4pew_1.4',
        # 'tip4pew_1.8', 'tip4pew_2.4', 'tip4pew_3']

        names = ["1", "1.4", "1.8", "2.4", "3"]
        group_number = 0
        # range is 15 because i have 5 groups with 3 models per tunnel, so 15 total boxes in the plot.
        for i in range(15):
            group_name = directories[index].split("_")[1] + "_" + str(names[group_number])
            # condition to check if the SCID is present in current comparative_analysis group
            if scid in scids_in_group[group_name][0][name_of_tunnel]:
                keys = directories[index:index + 5]
                combined = defaultdict(list)
                entry = defaultdict(list)
                release = defaultdict(list)
                for key in keys:
                    if not current_values[0][key]:  # if no value for the simulation, use 0
                        entry.setdefault(key,[])
                    else:
                        entry[key]= current_values[0][key]
                    if not current_values[1][key]:
                        release.setdefault(key,[])
                    else:
                        release[key]=current_values[1][key]
                    if combined_values_all[key] != 0:
                        combined[key] = combined_values_all[key]
                    else:
                        combined.setdefault(key,[])

                def update_dict(old_dict:dict, new_dict:dict):
                    for key,value in new_dict.items():
                        if key in old_dict:
                            old_dict[key].extend(new_dict[key])
                        else:
                            old_dict[key] = new_dict[key]

                update_dict(combined_dict,combined)
                update_dict(entry_dict,entry)
                update_dict(release_dict,release)

            else:
                if debug:
                    print(f"SCID - {scid} not present in {group_name}, not processing data from it")
                    keys = directories[index:index + 5]
                    for key in keys:
                        print(combined_values_all[key])
                else:
                    pass
            index += 5  # move to next 5 simulations

            # Move to next group / there are 15 simulations per water model and // 5 groups, which is why range(0,5)
            if index % 15 in list(range(0, 5)):
                group_number += 1

    return combined_dict, entry_dict, release_dict
