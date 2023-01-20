#  Find number of traced waters across
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

models = ['opc', 'tip3p', 'tip4pew']
epoch_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
sim_list = ["1", "2", "3", "4", "5"]
traced_waters = []
OPC_caver = (
(33242, 17928, 18804, 39839, 26756), (22777, 28839, 31777, 33444, 29527), (24562, 34476, 34924, 32485, 33225),
(41495, 34153
 , 28703, 33566, 29211), (77188, 74622, 29496, 52460, 45273))
TIP3P_caver = (
(33696, 19908, 20447, 40690, 28888), (21994, 28241, 27567, 34304, 26321), (25272, 34093, 35367, 28194, 32639), (45633,
                                                                                                                31160,
                                                                                                                26253,
                                                                                                                32041,
                                                                                                                37869),
(81293, 61061, 29398, 42491, 36490))
TIP4Pew_caver = (
(34560, 18193, 21899, 40603, 27051), (23532, 26750, 30350, 35084, 29328), (25620, 36190, 34411, 30850, 33471), (40686,
                                                                                                                34452,
                                                                                                                28903,
                                                                                                                32146,
                                                                                                                28344),
(71233, 70396, 28603, 40044, 38423))


def read_input():
    model_wat = []
    for model in models:
        epoch_wat = []
        for epoch in epoch_list:
            sim_wat = []
            for sim in sim_list:
                aq_results_dir = '/mnt/NAS/dean_water_models/md/' + model + '/' + epoch + '/' + sim + '/aquaduct/'
                with open(aq_results_dir + "5_analysis_results.txt") as txt:
                    read = txt.readlines()
                    # print(read)
                    txt.close()
                    line_number = read.index("Names of traced molecules: WAT\n")
                    sim_wat.append(int(read[line_number + 2].split(sep=':')[1]))
                    # traced_waters.append(model+','+epoch+','+sim+','+read[line_number+2].split(sep=':')[1])
            epoch_wat.append([*sim_wat])
        model_wat.append([*epoch_wat])
    return model_wat


def plot(water_data):
    width = 0.35
    import seaborn as sns
    sns.set_theme(style='ticks', font_scale=0.5)

    fig, ax = plt.subplots(1, 5, figsize=(12, 3), dpi=300)
    plt.subplots_adjust(wspace=0.2)
    fig.suptitle('Number of traced waters per simulation per epoch', fontsize=16, fontweight='bold')
    x = np.arange(1, 6)
    for epoch in range(5):
        ax[epoch].set_title(r'{}Å'.format(epoch_list[epoch][:-1]), size=14)
        # plt.show()
        ax[epoch].set_xticks(x)

        ax[epoch].yaxis.grid(True)
        opc_x = x - width / 3
        opc_height = water_data[0][epoch]
        tip3p_x = x
        tip3p_height = water_data[1][epoch]
        tip4pew_x = x + width / 3
        tip4pew_height = water_data[2][epoch]
        ax[epoch].bar(opc_x, opc_height, width / 3, label='OPC')
        ax[epoch].bar(tip3p_x, tip3p_height, width / 3, label='TIP3P')
        ax[epoch].bar(tip4pew_x, tip4pew_height, width / 3, label='TIP4P-Ew')
        ax[epoch].tick_params(axis='both', which='major', labelsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=14)
    ax[2].set_xlabel('# of simulation', fontsize=14)
    ax[0].set_ylabel('Number of water molecules', fontsize=14)
    plt.tight_layout()
    plt.savefig("/home/aravind/PhD_local/dean/figures/caver&aquaduct/aquaduct.png")
    # plt.show()


def plot_aquaduct_per_group(water_data):
    import seaborn as sns
    plt.style.use("grayscale")
    # sns.set(style='monochrome')
    # sns.color_palette('bright')
    fig, ax = plt.subplots(5, 3, figsize=(8.27, 11.7),dpi=300)
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
        return pd.DataFrame(dic,index=[1, 2, 3, 4, 5])
    opc_df = make_df(opc)
    tip3p_df = make_df(tip3p)
    tip4pew_df = make_df(tip4pew)
    # plt.subplots_adjust(wspace=0.4,hspace=0.3)
    plt.suptitle('AQUA-DUCT Traced Waters',fontsize=15,fontweight='bold')
    # fig.text(s='Sim Number',x = 0.48, y=0.07,fontsize=15,fontweight='bold')
    # fig.text(s='Number of waters', x=0.05, y=0.37, fontsize=15,rotation=90,fontweight='bold')
    opc_df.plot(ax=ax[0, 0], kind='bar', y='1A', legend=False, rot=0,
                color='r').set_ylabel("Group 1A",fontsize=15)
    opc_df.plot(ax=ax[1, 0], kind='bar', y='1.4A', legend=False, rot=0,
                color='r').set_ylabel("Group 1.4A",fontsize=15)
    opc_df.plot(ax=ax[2, 0], kind='bar', y='1.8A', legend=False,
                color='r').set_ylabel("Group 1.8A",fontsize=15)
    opc_df.plot(ax=ax[3, 0], kind='bar', y='2.4A', legend=False, rot=0,
                color='r').set_ylabel("Group 2.4A",fontsize=15)
    opc_df.plot(ax=ax[4, 0], kind='bar', y='3A', legend=False, rot=0,
                color='r').set_ylabel("Group 3A",fontsize=15)
    tip3p_df.plot(ax=ax[0, 1], kind='bar', y='1A', legend=False, rot=0,color='g')
    tip3p_df.plot(ax=ax[1, 1], kind='bar', y='1.4A', legend=False, rot=0,color='g')
    tip3p_df.plot(ax=ax[2, 1], kind='bar', y='1.8A', legend=False, rot=0,color='g')
    tip3p_df.plot(ax=ax[3, 1], kind='bar', y='2.4A', legend=False, rot=0,color='g')
    tip3p_df.plot(ax=ax[4, 1], kind='bar', y='3A', legend=False, rot=0,color='g')
    tip4pew_df.plot(ax=ax[0, 2], kind='bar', y='1A', legend=False, rot=0,color='c')
    tip4pew_df.plot(ax=ax[1, 2], kind='bar', y='1.4A', legend=False, rot=0,color='c')
    tip4pew_df.plot(ax=ax[2, 2], kind='bar', y='1.8A', legend=False, rot=0,color='c')
    tip4pew_df.plot(ax=ax[3, 2], kind='bar', y='2.4A', legend=False, rot=0,color='c')
    tip4pew_df.plot(ax=ax[4, 2], kind='bar', y='3A', legend=False, rot=0,color='c')
    ax[0, 0].set_title("P1", fontsize=15, fontweight='bold')
    ax[0, 1].set_title("P2", fontsize=15, fontweight='bold')
    ax[0, 2].set_title("P3", fontsize=15, fontweight='bold')
    ax[4,1].set_xlabel("Simulation Number",fontsize=10,fontweight='bold')
    for axes in ax.flatten():
        axes.tick_params(which='both',size=3,length=5)
    # plt.show()
    plt.tight_layout(pad=1, w_pad=0.5, h_pad=0.9)
    plt.savefig("/home/aravind/PhD_local/dean/figures/caver&aquaduct/waters_per_group.png")

    # sns.barplot(ax=ax[0,0],data=)


def plot_caver():
    width = 0.35
    import seaborn as sns
    sns.set_theme(style='white', font_scale=0.5)

    fig, ax = plt.subplots(1, 5, figsize=(12, 3), sharey='row', dpi=300)
    plt.subplots_adjust(wspace=0.2)
    fig.suptitle('Number of caver tunnels per simulation per epoch', fontsize=16, fontweight='bold')
    x = np.arange(1, 6)
    for epoch in range(5):
        ax[epoch].set_title(r'{}Å'.format(epoch_list[epoch][:-1]), size=14)
        # plt.show()
        ax[epoch].set_xticks(x)

        ax[epoch].yaxis.grid(True)
        opc_x = x - width / 3
        opc_height = OPC_caver[epoch]
        tip3p_x = x
        tip3p_height = TIP3P_caver[epoch]
        tip4pew_x = x + width / 3
        tip4pew_height = TIP4Pew_caver[epoch]
        ax[epoch].bar(opc_x, opc_height, width / 3, label='OPC')
        ax[epoch].bar(tip3p_x, tip3p_height, width / 3, label='TIP3P')
        ax[epoch].bar(tip4pew_x, tip4pew_height, width / 3, label='TIP4P-Ew')
        ax[epoch].tick_params(axis='both', which='major', labelsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=14)
    ax[2].set_xlabel('# of simulation', fontsize=14)
    ax[0].set_ylabel('Number of caver tunnels', fontsize=14)
    plt.tight_layout()
    plt.savefig("/home/aravind/PhD_local/dean/figures/caver.png")
    plt.show()


waters= read_input()
# plot_aquaduct_per_group(waters)
plot(waters)
# plot_caver()
