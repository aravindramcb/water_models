# -*- coding: utf-8 -*-
# Find if waters are forming H-bond networks within themselves and migrating as a cluster
import os
import pickle
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np
from matplotlib import pyplot as plt

@dataclass(order=True)
class TEvent:
    tid: int = field(compare=True, hash=True)
    watid: int = 0
    frame: int = -1
    radius: float = 999.999
    dist2lig: float = 999.999
    fraction: float = -1.0
    event_type: str = ""
    supercluster: int = -1
    caver_cluster: int = -1


def plot_radii_by_HB(data, hbond_outfolder, outfolder, tag="", frac_thr=0.7):
    def _get_hbonds(hbondfile):
        hbonds = 0
        with open(hbondfile, "r") as hfile:
            for line in hfile:
                if line.startswith("#"):
                    continue
                hbonds = int(line.strip().split()[1])  # 3 is event_H[Bridge]
        return hbonds

    mds = list(data.keys())
    mds.sort()
    for md in mds:
        for event in data[md]:
            # hbondfile = "{}/hbonds/{}/f{:0>5}_e{:0>4}.txt".format(outfolder, md, event.frame+1, event.tid)
            hbondfile = "{}/{}/f{:0>5}_e{:0>4}.txt".format(hbond_outfolder, md, event.frame + 1, event.tid)
            event.hbonds = _get_hbonds(hbondfile)

    sc_mds = defaultdict(int)
    for md in mds:
        _scs = set()
        for event in data[md]:
            if event.fraction >= frac_thr:
                _scs.add(event.supercluster)
        for s in _scs:
            sc_mds[s] += 1

    e_hb = defaultdict(list)
    for md in mds:
        for event in data[md]:
            # if water is buried in tunnel and SC present in more than one simulation
            if (event.fraction >= frac_thr) and (sc_mds[event.supercluster] >= 1):
                e_hb[event.hbonds].append(event.radius)

    ehb_hists = {}
    hbonds = list(e_hb.keys())
    hbonds.sort()
    for hb in hbonds:
        ehb_hists[hb], hb_xranges = np.histogram(e_hb[hb], range=(0.9, 3), bins=20)

    x_offset = (hb_xranges[1] - hb_xranges[0]) / 2.0

    bottoms = np.zeros((len(hbonds), ehb_hists[hbonds[0]].shape[0]))
    for i in range(1, len(hbonds)):
        for j in range(i, len(hbonds)):
            bottoms[j, :] += ehb_hists[hbonds[i - 1]]

    fig, ax = plt.subplots()
    for i in range(len(hbonds)):
        ax.bar(hb_xranges[:-1] + x_offset, height=ehb_hists[hbonds[i]], width=0.09,
               bottom=bottoms[i], label="H-bonds_{}".format(hbonds[i]))
    ax.set_xlim(hb_xranges[0] - 0.01, hb_xranges[-1] + 0.01)
    ax.set_ylim(0, 100)  # dhaa:300, epx:2250, lipase:1400
    ax.set_xticks(hb_xranges)
    ax.set_xticklabels(labels=[str(np.round(x, 1)) for x in hb_xranges],
                       fontdict={"fontsize": 10})
    ax.legend(loc="upper right", fontsize=8)
    plt.ylabel("Number of events")
    plt.xlabel(r"Radii [$\AA$]")
    plt.title("Radii of closest sphere to water transport by H-bonds in 1-MDs or more")
    plt.savefig(os.path.join(outfolder, tag + "radii_by_HB_bridge_md1.png"), dpi=300, format="png")


if __name__ == '__main__':
    dat_file = "/data/aravindramt/dean/tt/hbond/databases/database_1_8_T3.dat"
    with open(dat_file, "rb") as fin:
        database = pickle.load(fin)
    hbond_outfolder = "/data/aravindramt/dean/tt/hbond/hbonds/hb_1_8_T3"
    plot_folder = "/data/aravindramt/dean/tt/hbond/plots/1_8_TIP3P"

    plot_radii_by_HB(database, hbond_outfolder, plot_folder)
