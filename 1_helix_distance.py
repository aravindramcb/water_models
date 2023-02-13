# -*- coding: utf-8 -*-
# Measure and plot the CA distance between helices of tunnels.


import pandas as pd
import pytraj as pt
from matplotlib import pyplot as plt
import os

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


def collect_distance(results_dir, residues_to_load, atom_mask_to_measure):
    root_dir = results_dir
    distance_df = pd.DataFrame()
    for folder in sim_dir:
        trajectory = os.path.join(root_dir, folder, 'merged.nc')
        topology = os.path.join(root_dir, folder, 'structure_HMR.parm7')

        loaded_traj = pt.iterload(trajectory, topology, mask=residues_to_load)
        distance = pt.distance(loaded_traj, atom_mask_to_measure)
        print(distance)
        distance_df[folder] = distance
    return distance_df


# Plot openings
def plot(tunnel):
    import seaborn as sns
    distance_from_csv = pd.read_csv(f'/mnt/gpu/dean/md/{tunnel}_openings.csv', index_col='Unnamed: 0')
    distance_from_csv = distance_from_csv.reindex(columns=sim_dir)
    rows, cols = 5, 3
    fig, ax = plt.subplots(5, 3, figsize=(8.27, 11.7), sharex="col", sharey="row", dpi=300)
    plt.suptitle(f"{tunnel}".upper() + " Helix - Helix Distance (Å)", fontsize=15, fontweight='bold')
    i = 0
    j = 5
    x_label = ["Mean", 1, 2, 3, 4, 5]
    for row in range(rows):
        for col in range(cols):
            data = distance_from_csv.iloc[:, i:j]
            mean = {"Mean": data.mean(axis=1).values}
            mean_df = pd.DataFrame(data=mean)
            data = pd.concat([mean_df, data], axis=1)
            _plt = sns.boxplot(ax=ax[row, col], data=data, fliersize=0.2, width=0.6,
                               linewidth=0.5)
            _plt.artists[0].set_facecolor('grey')
            _plt.set_ylim(6, 16)
            # _plt.set_ylim(1, 19)
            # _plt.set_ylim(5, 10)
            ax[row, col].set_xticklabels(x_label)
            i += 5
            j += 5
    plt.tight_layout(pad=2.5, w_pad=0.5, h_pad=0.5)
    ax[0, 0].set_title("OPC", fontsize=15, fontweight='bold')
    ax[0, 1].set_title("TIP3P", fontsize=15, fontweight='bold')
    ax[0, 2].set_title("TIP4P-Ew", fontsize=15, fontweight='bold')
    ax[0, 0].set_ylabel("Group 1A", fontsize=15)
    ax[1, 0].set_ylabel("Group 1.4A", fontsize=15)
    ax[2, 0].set_ylabel("Group 1.8A", fontsize=15)
    ax[3, 0].set_ylabel("Group 2.4A", fontsize=15)
    ax[4, 0].set_ylabel("Group 3A", fontsize=15)
    plt.savefig(f'/home/aravind/PhD_local/dean/figures/helix-helix_distance/opening_{tunnel}_ind.png')


if __name__ == '__main__':
    p1mask = ':141@CA :172@CA'
    p2mask = ':242@CA :139@CA'
    p3mask = ':203@CA :149@CA'
    results = '/mnt/gpu/dean/md/simulations/'

    # MEASURE DISTANCE
    # distance_df = collect_distance(results_dir,residues_to_load,atom_mask_to_measure)
    distance_df = collect_distance(results, ':242,139', p2mask)
    # Always save as {your_name}_openings.csv
    distance_df.to_csv('./dean/md/p2_openings.csv')

    # PLOT SUBPLOTS OF MEASUREMENT
    plot('p1')
