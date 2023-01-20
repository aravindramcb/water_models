# Measure the CA distance between PHE141 and CYS171 (Main tunnel last residue)
import numpy as np
import pandas as pd
import pytraj as pt
from matplotlib import pyplot as plt
import os
import seaborn as sns
water_models = ['opc', 'tip3p', 'tip4pew']  #
epoch_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
sim_list = ["1", "2", "3", "4", "5"]
simulations = ['prod2.nc', 'prod3.nc']
topology = 'temp/structure_HMR.parm7'


def plot_average_bottleneck(model, epoch, simulation, distance_array):
    plt.style.use('seaborn-dark')
    fig, ax = plt.subplots()
    ax.margins(0.07)
    ax.set_title(f" {model} {epoch} {simulation}", fontsize=14)
    ax.set_xlabel("Frame Number", fontsize=13)
    ax.set_ylabel("Distance ($\AA$)", fontsize=13)
    # plt.grid(which='major', linestyle='-', color='darkorange')
    # plt.grid(True, which='minor', axis='both', linestyle='dashed', alpha=0.4, fillstyle='left', color='bisque')
    ax.plot(distance_array, label='Distance ($\AA$)')
    ax.legend()
    plt.show()
    print("plotting and saving")
    # plt.savefig("/home/aravind/PhD_local/dean/distance_main_tunnel/"+model+epoch+simulation+".png")
    # plt.close()


def multiplot(distance_of_all_systems):
    # plt.style.use('classic')
    colour_list = [["blue", "orange", "green", "red", "purple"],
                   ["purple", "red", "lawngreen", "orange", "blue"],
                   ["aqua", "fuchsia", "dodgerblue", "pink", "coral"]]
    y_tick = list(range(12, 19))
    rows, cols = 5, 5
    fig, ax = plt.subplots(5, 5, sharex='col', figsize=(25, 25))
    # plt.subplots_adjust(hspace=0.3, wspace=0.3)
    for row in range(rows):
        for col in range(cols):
            for system in [0, 1, 2]:
                distance_array = distance_of_all_systems[system]
                ax[row, col].set_title(f"{epoch_list[col]} - {sim_list[row]}", fontsize=15, name='Lucida sans unicode')
                # ax.set_xlabel("Frame Number", fontsize=13)
                # ax.set_ylabel("Distance ($\AA$)", fontsize=13)
                ax[row, col].plot(distance_array[col][row], label='Distance ($\AA$)', color=colour_list[system][col])
                plt.sca(ax[row, col])
                plt.yticks(y_tick)

    # ax.legend()
    fig.text(0.5, 0.025, 'Frames', ha='center', fontsize=15, )
    fig.text(0.03, 0.5, 'Distance ($\AA$)', fontsize=15, rotation='vertical')
    fig.suptitle("Time evolution of distances between CA of PHE141 and LYS172", weight='bold', fontsize=20)
    plt.tight_layout(pad=9, w_pad=0, h_pad=4)
    # plt.savefig(fname="/home/aravind/172.png")
    # plt.show()


root_dir = '/mnt/gpu/dean/md/simulations/'
sim_dir = [
    '1A_opc_1', '1A_opc_2', '1A_opc_3', '1A_opc_4', '1A_opc_5', '1A_tip3p_1',
    '1A_tip3p_2', '1A_tip3p_3', '1A_tip3p_4', '1A_tip3p_5', '1A_tip4pew_1', '1A_tip4pew_2', '1A_tip4pew_3',
    '1A_tip4pew_4', '1A_tip4pew_5', '1.4A_opc_1', '1.4A_opc_2', '1.4A_opc_3', '1.4A_opc_4', '1.4A_opc_5',
    '1.4A_tip3p_1', '1.4A_tip3p_2', '1.4A_tip3p_3', '1.4A_tip3p_4', '1.4A_tip3p_5',
    '1.4A_tip4pew_1', '1.4A_tip4pew_2', '1.4A_tip4pew_3', '1.4A_tip4pew_4', '1.4A_tip4pew_5', '1.8A_opc_1',
    '1.8A_opc_2', '1.8A_opc_3', '1.8A_opc_4', '1.8A_opc_5',
    '1.8A_tip3p_1', '1.8A_tip3p_2', '1.8A_tip3p_3', '1.8A_tip3p_4', '1.8A_tip3p_5',
    '1.8A_tip4pew_1', '1.8A_tip4pew_2', '1.8A_tip4pew_3', '1.8A_tip4pew_4', '1.8A_tip4pew_5',
    '2.4A_opc_1', '2.4A_opc_2', '2.4A_opc_3', '2.4A_opc_4', '2.4A_opc_5', '2.4A_tip3p_1', '2.4A_tip3p_2',
    '2.4A_tip3p_3', '2.4A_tip3p_4', '2.4A_tip3p_5', '2.4A_tip4pew_1', '2.4A_tip4pew_2', '2.4A_tip4pew_3',
    '2.4A_tip4pew_4', '2.4A_tip4pew_5', '3A_opc_1', '3A_opc_2', '3A_opc_3', '3A_opc_4', '3A_opc_5', '3A_tip3p_1',
    '3A_tip3p_2', '3A_tip3p_3', '3A_tip3p_4', '3A_tip3p_5', '3A_tip4pew_1', '3A_tip4pew_2', '3A_tip4pew_3',
    '3A_tip4pew_4', '3A_tip4pew_5']
# sim_dir = sorted(dirs)

distance_df = pd.DataFrame()
rmsf_df = pd.DataFrame()

def collect_distance(resids, mask):
    for folder in sim_dir:
        trajectory = os.path.join(root_dir, folder, 'merged.nc')
        topology = os.path.join(root_dir, folder, 'structure_HMR.parm7')

        loaded_traj = pt.iterload(trajectory, topology, mask=resids)
        distance = pt.distance(loaded_traj, mask)
        print(distance)
        distance_df[folder] = distance
def collect_rmsf():
    for folder in sim_dir:
        trajectory = os.path.join(root_dir, folder, 'merged.nc')
        topology = os.path.join(root_dir, folder, 'structure_HMR.parm7')
        print(trajectory, topology)
        loaded_traj = pt.iterload(trajectory, topology, mask='@CA,C,N')
        crystal = pt.load('/mnt/gpu/dean/md/simulations/crystal.pdb', top='/mnt/gpu/dean/md/simulations/crystal.parm7',
                          mask='@CA,C,N')
        pt.superpose(loaded_traj, ref=crystal)
        rmsf = pt.rmsf(loaded_traj, '@CA,C,N')
        rmsf_df[folder] = rmsf[:,1]


# collect_rmsf()
# rmsf_df.to_csv('./dean/md/rmsf.csv')
p1mask = ':141@CA :172@CA'
p2mask = ':242@CA :139@CA'
p3mask = ':203@CA :149@CA'
# collect_distance(':242,139', p2mask)
# distance_df.to_csv('./dean/md/p2_openings.csv')

# dump_file=os.path.join(root_dir,folder,'CA_distance.dump')
# dist=np.load(dump_file,allow_pickle=True)

# Plot openings
def plot(tunnel):
    import seaborn as sns
    distance_from_csv = pd.read_csv(f'/mnt/gpu/dean/md/{tunnel}_openings.csv', index_col='Unnamed: 0')
    distance_from_csv = distance_from_csv.reindex(columns=sim_dir)
    rows, cols = 5, 3
    fig, ax = plt.subplots(5, 3, figsize=(8.27, 11.7), sharex="col",sharey="row",dpi=300)
    plt.suptitle(f"{tunnel}".upper()+" Helix - Helix Distance (Å)", fontsize=15, fontweight='bold')
    i = 0
    j = 5
    x_label = ["Mean",1, 2, 3, 4, 5]
    for row in range(rows):
        for col in range(cols):
            data=distance_from_csv.iloc[:, i:j]
            mean={"Mean":data.mean(axis=1).values}
            mean_df = pd.DataFrame(data=mean)
            data=pd.concat([mean_df,data],axis=1)
            _plt = sns.boxplot(ax=ax[row, col], data=data,fliersize=0.2,width=0.6,
                               linewidth=0.5)
            _plt.artists[0].set_facecolor('grey')
            _plt.set_ylim(6, 16)
            # _plt.set_ylim(1, 19)
            # _plt.set_ylim(5, 10)
            ax[row, col].set_xticklabels(x_label)
            i += 5
            j += 5
    plt.tight_layout(pad=2.5,w_pad=0.5,h_pad=0.5)
    ax[0, 0].set_title("OPC", fontsize=15, fontweight='bold')
    ax[0, 1].set_title("TIP3P", fontsize=15, fontweight='bold')
    ax[0, 2].set_title("TIP4P-Ew", fontsize=15, fontweight='bold')
    ax[0, 0].set_ylabel("Group 1A", fontsize=15)
    ax[1, 0].set_ylabel("Group 1.4A", fontsize=15)
    ax[2, 0].set_ylabel("Group 1.8A", fontsize=15)
    ax[3, 0].set_ylabel("Group 2.4A", fontsize=15)
    ax[4, 0].set_ylabel("Group 3A", fontsize=15)
    plt.savefig(f'/home/aravind/PhD_local/dean/figures/helix-helix_distance/opening_{tunnel}_ind.png')
    # plt.show()
plot('p1')


def plot_rmsf(input_csv=str()):
    rmsf = pd.read_csv(input_csv,index_col='Atom_Numbers')
    rmsf = rmsf.reindex(columns=sim_dir)
    fig,axes = plt.subplots(75)
    print(len(axes))
    pic=rmsf.plot(ax=axes,subplots=True,figsize=(20,15),sharex=True)

    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    plt.show()


# plot_rmsf('/mnt/gpu/dean/md/rmsf.csv')