# get transport tools events per model

import os

import pandas
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Overall Results
# Create Filter
root_dir = "/run/user/1000/gvfs/sftp:host=gpu/data/aravindramt/dean/tt/tt_full_1_5_clustering"
df = pd.read_csv(os.path.join(root_dir + '/statistics/4-filtered_events_statistics.txt'),
                 index_col=False, skiprows=20, skipfooter=1, engine='python', skipinitialspace=True, na_values="-")
df = df.fillna(0)
# print('The keys are as follows:')
# display(df.keys())
# Define tunnels
tunnels_dict = {"P1": (1, 2, 4, 10, 11, 12), "P2": (3, 5, 6, 9, 17, 25, 27), "P3": (7, 8, 24)}
# Fix dtype
dfx = df.astype(
    {"SC_ID": int, "No_Sims": int, "Total_No_Frames": int, "Num_Events": int, "Num_entries": int, "Num_releases": int})

# Initial filter - by SC_ID
P1 = dfx.loc[(dfx["SC_ID"].isin(tunnels_dict["P1"])) & (dfx["No_Sims"] > 2) & (dfx["Num_Events"] > 1)]
P2 = dfx.loc[(dfx["SC_ID"].isin(tunnels_dict["P2"])) & (dfx["No_Sims"] > 2) & (dfx["Num_Events"] > 1)]
P3 = dfx.loc[(dfx["SC_ID"].isin(tunnels_dict["P3"])) & (dfx["No_Sims"] > 2) & (dfx["Num_Events"] > 1)]

# print("Original Data")
# display(P1,P2,P3)
# Final filter - by Num_sims, Num_frames and Num_events

P1_filter = {'SC_ID': P1["SC_ID"], "5_sim_eq": (P1["No_Sims"] / 75) * 5, "Avg_No_Frames": P1["Avg_No_Frames"],
             "Fraction_Sims": P1["No_Sims"] / 75, "Fraction_Frames": P1["Total_No_Frames"] / (75 * 2000)
    , "Fraction_events": P1["Num_Events"] / 287874}
P1_filter_df = pd.DataFrame(P1_filter)

P2_filter = {'SC_ID': P2["SC_ID"], "5_sim_eq": (P2["No_Sims"] / 75) * 5, "Avg_No_Frames": P2["Avg_No_Frames"],
             "Fraction_Sims": P2["No_Sims"] / 75, "Fraction_Frames": P2["Total_No_Frames"] / (75 * 2000)
    , "Fraction_events": P2["Num_Events"] / 287874}
P2_filter_df = pd.DataFrame(P2_filter)

P3_filter = {'SC_ID': P3["SC_ID"], "5_sim_eq": (P3["No_Sims"] / 75) * 5, "Avg_No_Frames": P3["Avg_No_Frames"],
             "Fraction_Sims": P3["No_Sims"] / 75, "Fraction_Frames": P3["Total_No_Frames"] / (75 * 2000)
    , "Fraction_events": P3["Num_Events"] / 287874}
P3_filter_df = pd.DataFrame(P3_filter)

# print("Proposed Filter")
master_filter = P1_filter_df.append(P2_filter_df.append(P3_filter_df)).sort_index()
# display(master_filter)


# P2["SC_ID","No_Sims","Num_Events"]
# P3["SC_ID","No_Sims","Num_Events"]
root_dir = "/run/user/1000/gvfs/sftp:host=gpu/data/aravindramt/dean/tt/tt_full_1_5_clustering/statistics/comparative_analysis/"
folders = [f for f in os.listdir(root_dir) if os.path.isdir(f)]
groups_df = {"1A": [], "1.4A": [], "1.8A": [], "2.4A": [], "3A": []}
model_names = ["opc", "tip3p", "tip4pew"]
# GROUP RESULTS BY groups_df
# ---------------------------
for group in groups_df.keys():
    #     print("Results of group by group name for " ,group)
    #     print("---------------------------------------------------------------------------------------------------------------")

    # Read results file and append the results as dataframes in a groups_df dictionary
    for model in model_names:
        df = pd.read_csv(os.path.join(root_dir + model + "_" + group[:-1] + '/4-filtered_events_statistics.txt'),
                         index_col=False, skiprows=20, skipfooter=1, engine='python', skipinitialspace=True,
                         na_values="-")
        df = df.fillna(0)
        sim_name = model + "_" + group
        # Uncomment to print raw results for all groups
        #         print(f"Filtered transport events for {sim_name}")
        groups_df[group].append(df)

# STAGE 1 FILTER BY SC_ID AND GROUP BY P1,P2 AND P3
# ----------------------------------------------------
stage_1_filtered = {"1A": [], "1.4A": [], "1.8A": [], "2.4A": [], "3A": []}
for group in groups_df:
    #     print(f"Results of tunnel_groups for {group}")
    #     print("--------------------------------------")

    for sim_df in groups_df[group]:
        #         print(sim_df.head())
        tunnels_dict = {"P1": (1, 2, 4, 10, 11, 12), "P2": (3, 5, 6, 9, 17, 25, 27), "P3": (7, 8, 24)}
        # Fix dtype
        dfx = sim_df.astype(
            {"SC_ID": int, "No_Sims": int, "Total_No_Frames": int, "Num_Events": int, "Num_entries": int,
             "Num_releases": int})
        total_events = dfx["Num_Events"].sum()

        # Filter - by SC_ID + Num_Events and group by tunnels
        P1 = dfx.loc[(dfx["SC_ID"].isin(tunnels_dict["P1"])) & (dfx["Num_Events"] >= 1)]
        P2 = dfx.loc[(dfx["SC_ID"].isin(tunnels_dict["P2"])) & (dfx["Num_Events"] >= 1)]
        P3 = dfx.loc[(dfx["SC_ID"].isin(tunnels_dict["P3"])) & (dfx["Num_Events"] >= 1)]

        # TO DISCUSS
        # ----------
        # Group - by Num_sims =5, Num_frames=5*2000 and Fraction_events=Num_events/Total_events

        P1_filter = {'SC_ID': P1["SC_ID"], "No_Sims": P1["No_Sims"], "Total_No_Frames": P1["Total_No_Frames"],
                     "Fraction_Sims": P1["No_Sims"] / 5, "Fraction_Frames": P1["Total_No_Frames"] / (5 * 20000)
            , "Fraction_events": P1["Num_Events"] / total_events}
        P1_filter_df = pd.DataFrame(P1_filter)

        P2_filter = {'SC_ID': P2["SC_ID"], "No_Sims": P2["No_Sims"], "Total_No_Frames": P2["Total_No_Frames"],
                     "Fraction_Sims": P2["No_Sims"] / 5, "Fraction_Frames": P2["Total_No_Frames"] / (5 * 20000)
            , "Fraction_events": P2["Num_Events"] / total_events}
        P2_filter_df = pd.DataFrame(P2_filter)

        P3_filter = {'SC_ID': P3["SC_ID"], "No_Sims": P3["No_Sims"], "Total_No_Frames": P3["Total_No_Frames"],
                     "Fraction_Sims": P3["No_Sims"] / 5, "Fraction_Frames": P3["Total_No_Frames"] / (5 * 20000)
            , "Fraction_events": P3["Num_Events"] / total_events}
        P3_filter_df = pd.DataFrame(P3_filter)
        sim_filtered = P1_filter_df.append(P2_filter_df.append(P3_filter_df)).sort_index()
        #         display(sim_filtered)
        stage_1_filtered[group].append(sim_filtered)

# Checked and OK

# STAGE 2 - FILTER EACH GROUP
# ----------------------------
P1_scs = [1, 2, 4, 10, 11, 12]
P2_scs = [3, 5, 6, 9, 17, 25, 27]
P3_scs = [7, 8, 24]
final_selection = []
for f1_group in stage_1_filtered:  # 1A,1.4A...
    #     print("Group", f1_group)
    #     print("----------------")
    _tunnels_selected = {"P1": [], "P2": [], "P3": []}
    _models = ("OPC", "TIP3P", "TIP4P-Ew")
    n = 0
    for model in stage_1_filtered[f1_group]:  # 1Aopc,1Atip3p...
        # Uncomment to display unfiltered fraction for P1,P2,P3
        # print(_models)

        # Filter master_filter to display only SC_ID in the current sim group (eg: opc_1.4A)
        filtered_master_df = master_filter[master_filter.index.isin(model.index)]

        # TO DISCUSS
        # ---------------
        # Filter by sims, frames
        # Current filter: Number of frames >= 1000 ,

        model["Fil_No_sims"] = np.where((model["No_Sims"] >= 1), "yes", "no")
        model["fil_No_frames"] = np.where((model["Total_No_Frames"] >= 2000), "yes", "no")
        # uncomment to activate events filter
        #         model["Fil_avg_frames"] = np.where((model["Avg_No_Frames"] >= filtered_master_df["Avg_No_Frames"])
        #                                           ,"yes","no")
        #         model["Fr_sim>5"] = np.where((model["Fraction_Sims"] > 0.5),"yes","no")
        #         model['Fil_Fr_events'] = np.where((model["Fraction_events"] >= filtered_master_df["Fraction_events"]),"yes","no")
        selected = model.loc[(model["Fil_No_sims"] == "yes") & (model["fil_No_frames"] == "yes")]
        #         & (model["Fil_avg_frames"] == "yes") |(model["Fr_sim>5"] == "yes")& (model["Fil_Fr_events"] == "yes"))]

        #         @@@@@ DEBUG @@@@@@@
        #         uncomment to print unfiltered df
        #         print("For "+_models[n])
        #         display(model)
        #         display(selected.iloc[:,:5])
        #         display("Master_filter values")
        #         display(filtered_master_df)
        #         print("")
        n += 1

        # Display results per tunnel group
        _tmp_selected_id = selected["SC_ID"].tolist()
        _tunnels_selected["P1"].append([x for x in _tmp_selected_id if x in P1_scs])
        _tunnels_selected["P2"].append([x for x in _tmp_selected_id if x in P2_scs])
        _tunnels_selected["P3"].append([x for x in _tmp_selected_id if x in P3_scs])
    #     print('Group',f1_group)
    #     display(pd.DataFrame(_tunnels_selected, index=['OPC','TIP3P','TIP4P-Ew']))
    #     print("")
    final_selection.append(_tunnels_selected)
# for group in final_selection:
#     display(pd.DataFrame(group))

# STAGE 3- STUDY STATISTICS WITH FILTERED GROUPS:
# ------------------------------------------------
# Super cluster IDs of tunnels after filtering
# opc_scid = [{"P1": [1], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [1], "P2": [3, 5], "P3": [8]},
#             {"P1": [1], "P2": [3, 5], "P3": [8, 24]}
#     , {"P1": [1, 2, 4, 10, 11], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [2, 4, 10, 12], "P2": [3, 9, 17, 27], "P3": []}]
#
# tip3p_scid = [{"P1": [1], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [1], "P2": [3, 5], "P3": [7, 8]},
#               {"P1": [1, 10], "P2": [3, 5], "P3": [7, 8]}
#     , {"P1": [1, 2, 4, 10, 11], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [2, 4, 10, 12], "P2": [3, 9, 27], "P3": []}]
#
# tip4pew_scid = [{"P1": [1], "P2": [3, 5, 6], "P3": [7, 8]}, {"P1": [1], "P2": [3, 5], "P3": [7, 8]},
#                 {"P1": [1], "P2": [3, 5], "P3": [8, 24]}
#     , {"P1": [1, 2, 10], "P2": [6, 9], "P3": [7, 8]}, {"P1": [2, 4, 10], "P2": [3, 9], "P3": []}]

opc_scid = [{"P1": [1,11], "P2": [3, 5, 6,25], "P3": [7, 8]}, {"P1": [1,10,11], "P2": [3, 5,25], "P3": [8]},
            {"P1": [1], "P2": [3, 5, 25], "P3": [7,8, 24]}
    , {"P1": [1, 2, 4, 10, 11,12], "P2": [3, 5, 6,17,25], "P3": [7, 8]},
            {"P1": [2, 4, 10, 12], "P2": [3,5,6, 9, 17,25,27], "P3": [8,24]}]

tip3p_scid = [{"P1": [1,10,11], "P2": [3, 5, 6,27], "P3": [7, 8,24]}, {"P1": [1,10,11], "P2": [3, 5,25], "P3": [7, 8]},
              {"P1": [1, 10,11], "P2": [3, 5,6,25,27], "P3": [7, 8, 24]}
    , {"P1": [1, 2, 4, 10, 11, 12], "P2": [3, 5, 6, 17, 25], "P3": [7, 8,24]},
              {"P1": [2, 4, 10, 11,12], "P2": [3,5,6, 9,17,25, 27], "P3": [7,8,24]}]

tip4pew_scid = [{"P1": [1,11], "P2": [3, 5, 6,25], "P3": [7, 8]},
                {"P1": [1,10], "P2": [3, 5,6,25], "P3": [7, 8]},
                {"P1": [1], "P2": [3, 5,25], "P3": [8, 24]},
                {"P1": [1, 2,4,10,11,12], "P2": [3,5,6, 9,17,25,27], "P3": [7, 8,24]},
                {"P1": [2, 4, 10,11,12], "P2": [3,5,6, 9,17,25,27], "P3": [7,8]}]

grp_names = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
mod_names = []
df_by_bottleneck = {}
opc_df = {"1A": [], "1.4A": [], "1.8A": [], "2.4A": [], "3A": []}
i = 0
keys = list(opc_df)

for group in opc_df:
    # iterating through every group in groups_df for OPC to get the value of P1,P2 and P3
    # <tunnel>_sc_ids = [sc_id]
    # print("For group", keys[i])
    opc = groups_df[grp_names[i]][0]  # [1] for OPC
    _p1_df = opc.loc[(opc["SC_ID"].isin(opc_scid[i]["P1"]))]
    _p2_df = opc.loc[(opc["SC_ID"].isin(opc_scid[i]["P2"]))]
    _p3_df = opc.loc[(opc["SC_ID"].isin(opc_scid[i]["P3"]))]
    num_events_p1 = _p1_df["Num_Events"].sum()
    num_events_p2 = _p2_df["Num_Events"].sum()
    num_events_p3 = _p3_df["Num_Events"].sum()
    opc_df[keys[i]] = [num_events_p1, num_events_p2, num_events_p3]
    i += 1
opc_num_events_df = pd.DataFrame(opc_df, index=['P1', 'P2', 'P3'])
i = 0
tip3p_df = {"1A": [], "1.4A": [], "1.8A": [], "2.4A": [], "3A": []}
for group in tip3p_df:
    # iterating through every group in groups_df for OPC to get the value of P1,P2 and P3
    # <tunnel>_sc_ids = [sc_id]
    # print("For group", keys[i])
    tip3p = groups_df[grp_names[i]][1]  # [1] for tip3p
    _p1_df = tip3p.loc[(tip3p["SC_ID"].isin(opc_scid[i]["P1"]))]
    _p2_df = tip3p.loc[(tip3p["SC_ID"].isin(opc_scid[i]["P2"]))]
    _p3_df = tip3p.loc[(tip3p["SC_ID"].isin(opc_scid[i]["P3"]))]
    num_events_p1 = _p1_df["Num_Events"].sum()
    num_events_p2 = _p2_df["Num_Events"].sum()
    num_events_p3 = _p3_df["Num_Events"].sum()
    tip3p_df[keys[i]] = [num_events_p1, num_events_p2, num_events_p3]
    i += 1
tip3p_num_events_df = pd.DataFrame(tip3p_df, index=['P1', 'P2', 'P3'])
i = 0
tip4pew_df = {"1A": [], "1.4A": [], "1.8A": [], "2.4A": [], "3A": []}
for group in tip4pew_df:
    # iterating through every group in groups_df for OPC to get the value of P1,P2 and P3
    # <tunnel>_sc_ids = [sc_id]
    # print("For group", keys[i])
    tip4pew = groups_df[grp_names[i]][2]  # [2] for tip4pew
    _p1_df = tip4pew.loc[(tip4pew["SC_ID"].isin(opc_scid[i]["P1"]))]
    _p2_df = tip4pew.loc[(tip4pew["SC_ID"].isin(opc_scid[i]["P2"]))]
    _p3_df = tip4pew.loc[(tip4pew["SC_ID"].isin(opc_scid[i]["P3"]))]
    num_events_p1 = _p1_df["Num_Events"].sum()
    num_events_p2 = _p2_df["Num_Events"].sum()
    num_events_p3 = _p3_df["Num_Events"].sum()
    tip4pew_df[keys[i]] = [num_events_p1, num_events_p2, num_events_p3]
    i += 1
tip4pew_num_events_df = pd.DataFrame(tip4pew_df, index=['P1', 'P2', 'P3'])

# TO PLOT CONSOLIDATED PLOT
# GROUP BY GROUP NAME AND TUNNEL
overall_df = [opc_num_events_df, tip3p_num_events_df, tip4pew_num_events_df]

# print(opc_num_events_df)

# GROUP BY MODEL AND TUNNEL
df_by_bottleneck['1A'] = pd.concat([opc_num_events_df.iloc[:, 0:1], tip3p_num_events_df.iloc[:, 0:1],
                                    tip4pew_num_events_df.iloc[:, 0:1]], axis=1)
df_by_bottleneck['1A'].columns = ['opc', 'tip3p', 'tip4pew']

df_by_bottleneck['1.4A'] = pd.concat([opc_num_events_df.iloc[:, 1:2], tip3p_num_events_df.iloc[:, 1:2],
                                      tip4pew_num_events_df.iloc[:, 1:2]], axis=1)
df_by_bottleneck['1.4A'].columns = ['opc', 'tip3p', 'tip4pew']

df_by_bottleneck['1.8A'] = pd.concat([opc_num_events_df.iloc[:, 2:3], tip3p_num_events_df.iloc[:, 2:3],
                                      tip4pew_num_events_df.iloc[:, 2:3]], axis=1)
df_by_bottleneck['1.8A'].columns = ['opc', 'tip3p', 'tip4pew']

df_by_bottleneck['2.4A'] = pd.concat([opc_num_events_df.iloc[:, 3:4], tip3p_num_events_df.iloc[:, 3:4],
                                      tip4pew_num_events_df.iloc[:, 3:4]], axis=1)
df_by_bottleneck['2.4A'].columns = ['opc', 'tip3p', 'tip4pew']

df_by_bottleneck['3A'] = pd.concat([opc_num_events_df.iloc[:, 4:], tip3p_num_events_df.iloc[:, 4:],
                                    tip4pew_num_events_df.iloc[:, 4:]], axis=1)
df_by_bottleneck['3A'].columns = ['opc', 'tip3p', 'tip4pew']


# df_by_bottleneck['1&1.4A'] = df_by_bottleneck['1A'].add(df_by_bottleneck['1.4A'], axis=1)
# print(df_by_bottleneck['1A'],df_by_bottleneck['1.4A'],df_by_bottleneck['1&1.4A'])

# GROUP BY MODEL AND TUNNEL AND SCALE USING TIP3P
#
def scale_by_tip3p():
    df_by_bottleneck['1A'] = pd.concat([opc_num_events_df.iloc[:, 0:1].divide(tip3p_num_events_df.iloc[:, 0:1])
                                           , tip3p_num_events_df.iloc[:, 0:1].divide(tip3p_num_events_df.iloc[:, 0:1]),
                                        tip4pew_num_events_df.iloc[:, 0:1].divide(tip3p_num_events_df.iloc[:, 0:1])],
                                       axis=1)
    df_by_bottleneck['1A'].columns = ['opc', 'tip3p', 'tip4pew']

    df_by_bottleneck['1.4A'] = pd.concat([opc_num_events_df.iloc[:, 1:2].divide(tip3p_num_events_df.iloc[:, 1:2]),
                                          tip3p_num_events_df.iloc[:, 1:2].divide(tip3p_num_events_df.iloc[:, 1:2]),
                                          tip4pew_num_events_df.iloc[:, 1:2].divide(tip3p_num_events_df.iloc[:, 1:2])],
                                         axis=1)
    df_by_bottleneck['1.4A'].columns = ['opc', 'tip3p', 'tip4pew']

    df_by_bottleneck['1.8A'] = pd.concat([opc_num_events_df.iloc[:, 2:3].divide(tip3p_num_events_df.iloc[:, 2:3]),
                                          tip3p_num_events_df.iloc[:, 2:3].divide(tip3p_num_events_df.iloc[:, 2:3]),
                                          tip4pew_num_events_df.iloc[:, 2:3].divide(tip3p_num_events_df.iloc[:, 2:3])],
                                         axis=1)
    df_by_bottleneck['1.8A'].columns = ['opc', 'tip3p', 'tip4pew']

    df_by_bottleneck['2.4A'] = pd.concat([opc_num_events_df.iloc[:, 3:4].divide(tip3p_num_events_df.iloc[:, 3:4]),
                                          tip3p_num_events_df.iloc[:, 3:4].divide(tip3p_num_events_df.iloc[:, 3:4]),
                                          tip4pew_num_events_df.iloc[:, 3:4].divide(tip3p_num_events_df.iloc[:, 3:4])],
                                         axis=1)
    df_by_bottleneck['2.4A'].columns = ['opc', 'tip3p', 'tip4pew']

    df_by_bottleneck['3A'] = pd.concat([opc_num_events_df.iloc[:, 4:].divide(tip3p_num_events_df.iloc[:, 4:]),
                                        tip3p_num_events_df.iloc[:, 4:].divide(tip3p_num_events_df.iloc[:, 4:]),
                                        tip4pew_num_events_df.iloc[:, 4:].divide(tip3p_num_events_df.iloc[:, 4:])],
                                       axis=1)
    df_by_bottleneck['3A'].columns = ['opc', 'tip3p', 'tip4pew']


# Plot final result without normalization
def plot_without_normalization():
    import matplotlib
    sns.set(style='white')
    matplotlib.rc('xtick', labelsize=10)
    matplotlib.rc('ytick', labelsize=12)
    matplotlib.rc('axes', titlesize=15)
    fig, ax = plt.subplots(5, 3, figsize=(8.27, 11.7),sharex='all')
    plt.suptitle("Events - Consolidated",fontsize=15,fontweight='bold')
    print(df_by_bottleneck["1A"].iloc[0:1, :])
    # P1
    g1a_plot = sns.barplot(ax=ax[0, 0], data=df_by_bottleneck["1A"].iloc[0:1, :],errorbar="sd").set(ylabel='GROUP 1A')
    g14_plot = sns.barplot(ax=ax[1, 0], data=df_by_bottleneck["1.4A"].iloc[0:1, :]).set(ylabel='GROUP 1.4A')
    # g1_plot=sns.boxplot(ax=ax[1,0],data=df_by_bottleneck["1&1.4A"].iloc[0:1,:]).set(ylabel='GROUP 1&1.4A')
    g18_plot = sns.barplot(ax=ax[2, 0], data=df_by_bottleneck["1.8A"].iloc[0:1, :]).set(ylabel='GROUP 1.8A')
    g24_plot = sns.barplot(ax=ax[3, 0], data=df_by_bottleneck["2.4A"].iloc[0:1, :]).set(ylabel='GROUP 2.4A')
    g3_plot = sns.barplot(ax=ax[4, 0], data=df_by_bottleneck["3A"].iloc[0:1, :]).set(ylabel='GROUP 3A')
    # P2
    g1a1_plot = sns.barplot(ax=ax[0, 1], data=df_by_bottleneck["1A"].iloc[1:2, :])
    g141_plot = sns.barplot(ax=ax[1, 1], data=df_by_bottleneck["1.4A"].iloc[1:2, :])
    # g11_plot=sns.barplot(ax=ax[1,1],data=df_by_bottleneck["1&1.4A"].iloc[1:2,:]).set(ylabel='GROUP 1&1.4A',title='. ')
    g181_plot = sns.barplot(ax=ax[2, 1], data=df_by_bottleneck["1.8A"].iloc[1:2, :])
    g241_plot = sns.barplot(ax=ax[3, 1], data=df_by_bottleneck["2.4A"].iloc[1:2, :])
    g31_plot = sns.barplot(ax=ax[4, 1], data=df_by_bottleneck["3A"].iloc[1:2, :])
    # P3
    g1a2plot = sns.barplot(ax=ax[0, 2], data=df_by_bottleneck["1A"].iloc[2:, :])
    g142_plot = sns.barplot(ax=ax[1, 2], data=df_by_bottleneck["1.4A"].iloc[2:3, :])
    # g112_plot=sns.barplot(ax=ax[1,2],data=df_by_bottleneck["1&1.4A"].iloc[2:3,:]).set(ylabel='GROUP 1&1.4A',title='. ')
    g182_plot = sns.barplot(ax=ax[2, 2], data=df_by_bottleneck["1.8A"].iloc[2:3, :])
    g242_plot = sns.barplot(ax=ax[3, 2], data=df_by_bottleneck["2.4A"].iloc[2:3, :])
    g32_plot = sns.barplot(ax=ax[4, 2], data=df_by_bottleneck["3A"].iloc[2:3, :])

    for axes in ax.flatten():
        l = axes.get_ylabel()
        x = axes.get_xlabel()
        axes.set_ylabel(l, fontsize=12, fontweight='bold')
        axes.set_xlabel(x, fontsize=12, fontweight='bold')
        axes.set_xticklabels(labels=['OPC','TIP3P','TIP4P-Ew'])
    # plt.subplots_adjust(left=0.1, right=0.902, bottom=0.056, top=0.92, wspace=0.202, hspace=0.3)
    plt.tight_layout(pad=1.5, w_pad=0.0, h_pad=0.2)
    h,l = ax[4,1].get_legend_handles_labels()
    ax[4,1].legend(h, ["OPC","TIP3P","TIP4P-Ew"])
    ax[4,2].set_yticklabels([])
    ax[0, 0].set_title("P1", fontweight='bold')
    ax[0, 1].set_title("P2", fontweight='bold')
    ax[0,2].set_title("P3",fontweight='bold')
    plt.savefig("/home/aravind/PhD_local/dean/figures/transport_tools/events.png")
    # plt.show()
    plt.close()



# COMBINED PLOTS
def normalized_plots():
    # NORMALIZE DATA AND PLOT - number of transport events
    # P1, P2, P3 before normalization
    p1_df_full = pd.concat([overall_df[0].iloc[0:1, :], overall_df[1].iloc[0:1, :], overall_df[2].iloc[0:1, :]], axis=0,
                           keys=["OPC", "TIP3P", "TIP4P-Ew"]).reset_index(level=1, drop=True)
    p2_df_full = pd.concat([overall_df[0].iloc[1:2, :], overall_df[1].iloc[1:2, :], overall_df[2].iloc[1:2, :]], axis=0,
                           keys=["OPC", "TIP3P", "TIP4P-Ew"]).reset_index(level=1, drop=True)
    p3_df_full = pd.concat([overall_df[0].iloc[2:, :], overall_df[1].iloc[2:, :], overall_df[2].iloc[2:, :]], axis=0,
                           keys=["OPC", "TIP3P", "TIP4P-Ew"]).reset_index(level=1, drop=True)
    opc_df_full = pd.concat([overall_df[0].iloc[:, 0:1], overall_df[1].iloc[:, 0:1], overall_df[2].iloc[:, 0:1]],
                            axis=0)
    sns.set(style='white', font_scale=1.3)
    fig, ax = plt.subplots(1, 3, figsize=(20, 5), dpi=300)
    plt.subplots_adjust(wspace=0.2)
    p1_plot = p1_df_full.plot(kind='bar', ax=ax[0], title="P1", rot=0).legend(loc='upper left')
    p2_plot = p2_df_full.plot(kind='bar', ax=ax[1], title="P2", rot=0)
    p3_plot = p3_df_full.plot(kind='bar', ax=ax[2], title="P3", rot=0)
    plt.tick_params(axis='both', which='major', labelsize=20)

    plt.savefig("/home/aravind/PhD_local/dean/figures/by_tunnel_before_normalization.png")
    plt.close()

    # https://datagy.io/pandas-normalize-column/
    def min_max_scaling(series):
        return (series - series.min()) / (series.max() - series.min())

    def z_score_standardization(series):
        return (series - series.mean()) / series.std()

    for dfs in [p1_df_full, p2_df_full, p3_df_full]:
        for col in dfs.columns:
            dfs[col] = z_score_standardization(dfs[col])

    print(p1_df_full)
    # P1, P2, P3 after normalization
    fig, ax = plt.subplots(1, 3, figsize=(20, 5), dpi=300)
    plt.subplots_adjust(top=0.85, wspace=0.2)
    plt.suptitle("Z Score Normalized")
    p1_plot = p1_df_full.plot(kind='bar', ax=ax[0], title="P1", rot=0, ylim=(-1, 1.5))
    p2_plot = p2_df_full.plot(kind='bar', ax=ax[1], title="P2", rot=0, ylim=(-1, 1.5))
    p3_plot = p3_df_full.plot(kind='bar', ax=ax[2], title="P3", rot=0, ylim=(-1, 1.5))
    plt.savefig("/home/aravind/PhD_local/dean/figures/by_tunnel_z_normalized.png")

if __name__ == '__main__':
    plot_without_normalization()