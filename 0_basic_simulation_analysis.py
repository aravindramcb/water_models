# basic simulation analyses script

import matplotlib.pyplot as plt
import pytraj as pt
import numpy as np
import pandas as pd
import os
import seaborn as sns
models = ['opc', 'tip3p', 'tip4pew']
epoch_list = ["1A", "1.4A", "1.8A", "2.4A", "3A"]
sim_list = ["1", "2", "3", "4", "5"]
# simulations = ['prod2.nc', 'prod3.nc']


def generate_rmsd(trajectory, topology):
    traj = pt.iterload(trajectory, topology, mask='!@WAT')
    rmsd = pt.rmsd(traj, mask=("@CA,@N,@C"), dtype='dataset')
    rmsf = pt.rmsf(traj, mask=("@CA,@N,@C"), options='byres')
    # Save the files
    df_rmsd = pd.DataFrame(rmsd)
    df_rmsf = pd.DataFrame(rmsf)
    df_rmsd.to_csv("./rmsd.dat", index=False)
    df_rmsf.to_csv("./rmsf.dat", index=False)
    return rmsd, rmsf


def process_rmsds(dat_file):
    with open(dat_file) as rmsd_data:
        data = rmsd_data.readlines()
        rmsd_data.close()
        str_rmsd = data[1].split(sep=',')
        rmsd_list = []
        time_list = []
        i = 0
        for line in str_rmsd:
            rmsd = float(line)
            rmsd_list.append(rmsd)
            time = float(i * 0.01)  # convert frame number to time in ns
            time_list.append(time)
            i += 1
    rmsd_array = np.array(rmsd_list).reshape(len(rmsd_list), 1)
    time_array = np.array(time_list)
    return rmsd_array, time_list
def process_rmsdf():
    with open('./rmsf.dat') as rmsf_data:
        data = rmsf_data.read()
        data_cleaned=data.splitlines()
        rmsf_data.close()
        rmsf_list = []
        time_list = []
        i = 0
        for line in data_cleaned:
            str_rmsf = line.split(sep=',')
            rmsf = float(str_rmsf[1])
            rmsf_list.append(rmsf)
            time = int(i)  # convert frame number to atom or residue
            time_list.append(time)
            i += 1
    rmsf_array = np.array(rmsf_list).reshape(len(rmsf_list), 1)
    time_array = np.array(time_list)

    return rmsf_array, time_list

def plot_rmsd(rmsd, time_scale):
    for sim in rmsd:
        # plt.figure(figsize=(5, 5))
        plt.plot(time_scale[1][1:], sim[1:])
        plt.xlabel('Time [ns]')
        plt.ylabel('RMSD [Å]')
    # plt.show()
        # plt.savefig('rmsd_backbone.png', dpi=300)

def plot_boxplot(data,systemname):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10,5), edgecolor="black", dpi=300)
    ax = sns.boxplot(data=data, fliersize=0.01, palette="Set2")
    # fig,ax=plt.plot(data)
    ax.set_ylim(0.3, 0.6)
    ax.grid(b=True, axis='x', linestyle='--')
    ax.set_title(systemname, fontsize=18)
    # ax.set_title(f"RMSF - {systemname}", fontsize=18)
    ax.set_xticklabels([str(i) for i in range(1, 6)], fontsize=10)
    ax.set_xlabel("Simulation Number", fontsize=16)
    ax.set_ylabel(f"RMSD ($\AA$)", fontsize=16)
    # ax.set_ylabel(f"RMSF ($\AA$)", fontsize=16)
    # plt.show()
    plt.savefig(f"/home/aravind/PhD_local/dean/figures/simulation_analysis/rmsd_{systemname}.png")
    plt.close()

def plot_subplot(data,bottleneck):
    sns.set()
    width=0.35
    fig,ax = plt.subplots(3,1,figsize=(8,10),sharex='row')
    plt.subplots_adjust(wspace=0.2)
    plt.suptitle(bottleneck,fontsize=15)
    ax[0].set_ylim(bottom=0.35, top=0.55)
    ax[1].set_ylim(bottom=0.35, top=0.55)
    ax[2].set_ylim(bottom=0.35,top=0.55)
    x = [1,2,3,4,5]
    a=sns.violinplot(ax=ax[0],data=data[0])
    b=sns.violinplot(ax=ax[1], data=data[1])
    p=sns.violinplot(ax=ax[2],data=data[2])
    a.set_ylabel("OPC", fontsize=15)
    b.set_ylabel("TIP3P", fontsize=15)
    p.set_ylabel("TIP4P-Ew", fontsize=15)
    p.set_xlabel("Simulation Number", fontsize=15)
    p.set_xticks(range(len(data[0])))
    p.set_xticklabels(x)
    a.set_xticklabels(x)
    b.set_xticklabels(x)
    plt.tight_layout()
    plt.show()
    # plt.savefig(f"/home/aravind/PhD_local/dean/figures/simulation_analysis/{bottleneck}_all_models.png")
def plot_violinplot(data,systemname):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 5), edgecolor="black", dpi=200)
    ax = sns.violinplot(data=data, palette="Set2")
    # fig,ax=plt.plot(data)
    ax.set_ylim(0.3, 0.6)
    ax.grid(b=True, axis='x', linestyle='--')
    ax.set_title(systemname, fontsize=18)
    # ax.set_title(f"RMSF - {systemname}", fontsize=18)
    ax.set_xticklabels([str(i) for i in range(1, 26)], fontsize=10)
    ax.set_xlabel("Simulation Number", fontsize=16)
    ax.set_ylabel(f"RMSD ($\AA$)", fontsize=16)
    # ax.set_ylabel(f"RMSF ($\AA$)", fontsize=16)
    # plt.show()
    plt.savefig(f"/home/aravind/PhD_local/dean/figures/simulation_analysis/rmsd_{systemname}_vp.png")
    plt.close()

def main():
    traj = './merged.nc'
    top_opc = 'temp/structure_HMR.parm7'
    top = 'structure_HMR.parm7'
    root_dir= '/mnt/gpu/dean/md/simulations'

    time_whole = []
    rmsf_whole=[]
    # 0=OPC, 1=TIP3p, 2=TIP4P-Ew
    os.chdir(root_dir)
    for epoch in epoch_list:
        rmsd_per_model= []
        for model in models:
            rmsd_per_bottleneck = []
            for simulation in sim_list:
                current_dir=epoch+"_"+model+"_"+simulation
                # os.chdir(current_dir)
                print('working in', current_dir)
                # generate_rmsd(traj,top)
                # rmsf_array, time = process_rmsdf()
                rmsd_array, time = process_rmsds(current_dir+'/rmsd.dat')
                # time_whole.append(time)
                # rmsf_whole.append(rmsf_array)
                rmsd_per_bottleneck.append(rmsd_array[1:])
            # plot_boxplot(rmsd_per_bottleneck,model+epoch)
            rmsd_per_model.append(rmsd_per_bottleneck)
        plot_subplot(rmsd_per_model,"RMSD-"+epoch)
    # for sim in rmsf_whole:
    #     plt.plot(time_whole[1][1:],sim[1:])
    #     plt.title('TIP3P')
    #     plt.xlabel('RES #')
    #     plt.ylabel('RMSf [Å2]')
    # plt.savefig("/home/aravind/PhD_local/dean/figures/simulation_analysis/rmsf_line_t3p.png",dpi=300)
    # plot_rmsd(rmsf_whole, time_whole)
    # plot_violinplot(rmsd_whole,'TIP4P-Ew')

if __name__ == '__main__':
    main()



