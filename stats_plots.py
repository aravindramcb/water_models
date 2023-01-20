# Load bottleneck values of SC <ID> from original caver calculations
# P1 tunnel
# [1,2,4,8,7,10]

# Load data from tunnel characteristics
import glob
import pandas as pd
import os
sim_dir = '/mnt/gpu/dean/tt/minimal_data/'
bottleneck_radius = pd.DataFrame()


def parse_bottlenecks(system, caver_cluster_id):
    tc = pd.read_csv(sim_dir + system + '/caver_analyses/final_clustering/analysis/tunnel_characteristics.csv',
                     usecols=[1, 5])
    tc_filtered = tc.loc[(tc[" Tunnel cluster"] == caver_cluster_id)]

    bottleneck_radius[system] = (tc_filtered[" Bottleneck radius"].reset_index(drop=True))


# sim_list = [i for i in os.listdir(sim_dir) if os.path.isdir(sim_dir + i)]
sim_list = ['1.4A_opc_1', '1.4A_opc_2', '1.4A_opc_3', '1.4A_opc_4', '1.4A_opc_5', '1.4A_tip3p_1', '1.4A_tip3p_2',
            '1.4A_tip3p_3', '1.4A_tip3p_4', '1.4A_tip3p_5', '1.4A_tip4pew_1', '1.4A_tip4pew_2', '1.4A_tip4pew_3',
            '1.4A_tip4pew_4', '1.4A_tip4pew_5', '1.8A_opc_1', '1.8A_opc_2', '1.8A_opc_3', '1.8A_opc_4', '1.8A_opc_5',
            '1.8A_tip3p_1', '1.8A_tip3p_2', '1.8A_tip3p_3', '1.8A_tip3p_4', '1.8A_tip3p_5', '1.8A_tip4pew_1',
            '1.8A_tip4pew_2', '1.8A_tip4pew_3', '1.8A_tip4pew_4', '1.8A_tip4pew_5', '1A_opc_1', '1A_opc_2', '1A_opc_3',
            '1A_opc_4', '1A_opc_5', '1A_tip3p_1', '1A_tip3p_2', '1A_tip3p_3', '1A_tip3p_4', '1A_tip3p_5', '1A_tip4pew_1',
            '1A_tip4pew_2', '1A_tip4pew_3', '1A_tip4pew_4', '1A_tip4pew_5', '2.4A_opc_1', '2.4A_opc_2', '2.4A_opc_3',
            '2.4A_opc_4', '2.4A_opc_5', '2.4A_tip3p_1', '2.4A_tip3p_2', '2.4A_tip3p_3', '2.4A_tip3p_4', '2.4A_tip3p_5',
            '2.4A_tip4pew_1', '2.4A_tip4pew_2', '2.4A_tip4pew_3', '2.4A_tip4pew_4', '2.4A_tip4pew_5', '3A_opc_1',
            '3A_opc_2', '3A_opc_3', '3A_opc_4', '3A_opc_5', '3A_tip3p_1', '3A_tip3p_2', '3A_tip3p_3', '3A_tip3p_4',
            '3A_tip3p_5', '3A_tip4pew_1', '3A_tip4pew_2', '3A_tip4pew_3', '3A_tip4pew_4', '3A_tip4pew_5']

# sorted_list = sorted(sim_list)
# print(sorted_list)
caver_id_p2= [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 3, 2, 2, 2, 3, 3, 2, 2, 2, 4, 3, 2, 7, 5, 2, 10, 2,
           3, 5, 2, 8, 2, 2, 6, 2, 6, 3, 2, 6, 9, 9, 3, 2, 2, 8, 9, 2, 10, 14, 16, 22, 2, 3, 2, 11, 6, 2, 2, 2, 7, 4, 2,
           2, 2, 25, 45]
# p3 - 7,8,24
caver_id_p3=[5, 4, 6, 5, 6, 3, 4, 5, 5, 5, 3, 4, 5, 5, 6, 4, 4, 8, 7, 6, 6, 4, 7, 6, 5, 5, 4, 7, 6, 6, 4, 4, 8, 8, 3,
            4,5, 9, 7, 3, 4, 4, 9, 6, 3, 0, 0, 10, 5, 4, 29, 0, 15, 4, 6, 35, 22, 10, 4, 5, 0, 0, 26, 0, 48, 0, 0, 18,
             31, 0, 0, 0, 17, 0, 0]
# caver_id_p3 = [5, 4, 6, 5, 6, 3, 4, 5, 5, 5, 3, 4, 5, 5, 6, 4, 4, 8, 7, 6, 6, 4, 7, 6, 5, 5, 4, 7, 6, 6, 4, 4, 8, 8, 3,
#                4, 5, 9, 7, 3, 4, 4, 9, 6, 3, 28, 14, 10, 5, 4, 29, 12, 8, 4, 6, 35, 22, 10, 4, 5, 19, 33, 26, 17, 48,
#                20, 23, 18, 31, 18, 18, 17, 17, 24, 29]

for sim in sim_list:
    index=list.index(sim_list,sim)
    print(sim,"->",caver_id_p3[index])
    parse_bottlenecks(sim,caver_id_p3[index])
    # parse_bottlenecks(sim,1)
bl_sorted_df = (bottleneck_radius.reindex(sorted(bottleneck_radius.columns), axis=1))

# Prepare group 1A and 1.4A
grp_1_4A = bl_sorted_df.iloc[:, 0:15]
grp_1_4A = grp_1_4A.dropna()
grp_1A = bl_sorted_df.iloc[:, 30:45]
grp_1A = grp_1A.dropna()
to_test = pd.concat([grp_1A, grp_1_4A], axis=1)

# tukey HSD test
# check documentation here https://github.com/reneshbedre/bioinfokit
# https://www.reneshbedre.com/blog/anova.html#one-way-one-factor-anova-with-python
# from bioinfokit.analys import stat
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# COMPARE 1A AND 1.4A
# ---------------------
# prepare data
# col_names=to_test.columns.values.tolist()
# df_melt=pd.melt(to_test.reset_index(),id_vars=['index'], value_vars=col_names)
# df_melt.columns=['index','Sim_ID','Bottleneck']
# fig,ax=plt.subplots(figsize=(35,10))
# plot=sns.boxplot(x='Sim_ID', y='Bottleneck', data=df_melt)
# plot.set_xticklabels(ax.get_xticklabels(),rotation=30)
# plt.show()

# COMPARE ALL TO ALL
# ------------------
# prepare data
col_names = bl_sorted_df.columns.values.tolist()
df_melt = pd.melt(bl_sorted_df.reset_index(), id_vars=['index'], value_vars=col_names)
df_melt.columns = ['index', 'Sim_ID', 'Bottleneck']

# Perform tukey_hsd test
# res = stat()
# res.tukey_hsd(df=df_melt, res_var='Bottleneck', xfac_var='Sim_ID', anova_model='Bottleneck ~ C(Sim_ID)')
# result = res.tukey_summary
# # result.to_excel("/mnt/gpu/dean/tt/tt_full_1_5_clustering/tukey_1_1_4.xlsx")
# result.to_excel("/mnt/gpu/dean/tt/tt_full_1_5_clustering/tukey_all_to_all.xlsx")

# Perform ANOVA using bioinfokit
# res=stat()
# res.anova_stat(df=df_melt, res_var="Bottleneck", anova_model='Bottleneck ~ C(Sim_ID)')
# res.anova_summary

# Prepare data for plotting
# Average bottleneck per group

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

# Concat data per group
df_1 = pd.concat([avg_1.loc[:, 'Whole'], bl_sorted_df.iloc[:, 30:45]], axis=1)
df_1_4 = pd.concat([avg_1_4.loc[:, 'Whole'], bl_sorted_df.iloc[:, 0:15]], axis=1)
df_1_8 = pd.concat([avg_1_8.loc[:, 'Whole'], bl_sorted_df.iloc[:, 15:30]], axis=1)
df_2_4 = pd.concat([avg_2_4.loc[:, 'Whole'], bl_sorted_df.iloc[:, 45:60]], axis=1)
df_3 = pd.concat([avg_3.loc[:, 'Whole'], bl_sorted_df.iloc[:, 60:75]], axis=1)

# Mean and Std
print("mean",avg_df.mean())
print("std",avg_df.std())
print("min", avg_df.min())
print("max",avg_df.max())

# ------------------------------------PLOT------------------------#


# Visualize Bottleneck as boxplot
import matplotlib
x_labels = ['Whole'] + [col_names[i].split(sep='_', maxsplit=1)[1] for i in range(15)]
fig, axes = plt.subplots(3, 2, figsize=(30, 25), dpi=150)
plt.subplots_adjust(hspace=0.3,top=1)
sns.set()
# plt.suptitle("BOTTLENECK RADII OF MAIN TUNNEL(P1)", fontsize=20, fontweight='bold')
plt.suptitle("BOTTLENECK RADII OF P3 TUNNEL", fontsize=20, fontweight='bold',y=0.98)
box1 = sns.boxplot(data=df_1, ax=axes[0, 0], color='b')
box1.set_xlabel("Group 1A", fontsize=20, fontweight='bold')
box1.artists[0].set_facecolor('grey')

box2 = sns.boxplot(data=df_1_4, ax=axes[0, 1], color='g')
box2.set_xlabel("Group 1.4A", fontsize=20, fontweight='bold')
box2.artists[0].set_facecolor('grey')

box3 = sns.boxplot(data=df_1_8, ax=axes[1, 0], color='r')
box3.set_xlabel("Group 1.8A", fontsize=20, fontweight='bold')
box3.set_ylabel("Bottleneck radii", fontsize=20, fontweight='bold')
box3.artists[0].set_facecolor('grey')

box4 = sns.boxplot(data=df_2_4, ax=axes[1, 1], color='c')
box4.set_xlabel("Group 2.4A", fontsize=20, fontweight='bold')
box4.artists[0].set_facecolor('grey')

box5=sns.boxplot(data=df_3,ax=axes[2,0],color='m')
box5.set_xlabel("Group 3A", fontsize=20, fontweight='bold')
box5.artists[0].set_facecolor('grey')
x_tick_location=np.arange(16)
for ax in axes.flatten():
    ax.set_ylim(bottom=0.8, top=3)
    ax.set_xticks(x_tick_location)
    ax.set_xticklabels(x_labels,rotation=30,fontsize=20)
    ax.tick_params(axis='y',labelsize=20)
    sns.set_theme(style='darkgrid')
color_pal = {'1A': 'b', '1.4A': 'g', '1.8A': 'r', '2.4A': 'c', '3A': 'm'}
box6 = sns.boxplot(data=avg_df,palette=color_pal)
box6.set_xlabel("Overall", fontsize=20, fontweight='bold')
plt.tight_layout(pad=1.8)
plt.savefig("/home/aravind/PhD_local/dean/figures/bottlenecks/p3.png")
# plt.show()
