# Python script to get the exact matching text file name for events having radii of 1A or lesser

import os
import pickle
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from collections import defaultdict
from dataclasses import dataclass, field

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

def get_narrow_tunnel_radii(data,hbond_outfolder, outfolder, tag="", frac_thr=0.7):
    def _get_hbonds(hbondfile):
        hbonds = 0
        with open(hbondfile, "r") as hfile:
            for line in hfile:
                if line.startswith("#"):
                    continue
                hbonds = int(line.strip().split()[2])
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


    file_names =defaultdict(list)
    for md in mds:
        for event in data[md]:
            if (event.fraction >= frac_thr) and (sc_mds[event.supercluster] >= 1) and (event.radius <= 1)\
                    and (event.hbonds >=1):
                hbondfile = "{}/{}/f{:0>5}_e{:0>4}.txt".format(hbond_outfolder, md, event.frame + 1, event.tid)
                exact_matching_file_name = f"{event.tid}_{event.event_type}_sc{event.supercluster}.txt"
                if md not in file_names:
                    file_names[md]=[[event.watid,hbondfile,exact_matching_file_name]]
                else:
                    file_names[md].append([event.watid,hbondfile,exact_matching_file_name])
    for key,value in file_names.items():
        if len(value) ==1:
            print(key, "-->", *value,"\n")
        else:
            for v in value:
                print(key, "-->", v, "\n")


if __name__ == '__main__':
    dat_file = "/data/aravindramt/dean/tt/hbond/databases/database_1_T3.dat"
    with open(dat_file, "rb") as fin:
        database = pickle.load(fin)
    hbond_outfolder = "/data/aravindramt/dean/tt/hbond/hbonds/hb_1_T3"
    plot_folder = "/data/aravindramt/dean/tt/hbond/plots/1_TIP3P"

    get_narrow_tunnel_radii(database, hbond_outfolder, plot_folder)