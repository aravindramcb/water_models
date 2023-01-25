#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import os
import pickle
import tempfile
import subprocess
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

def _ncontacts(mdtag, events, trajfolder, outfolder, prot_size):
    for event in events:
        ncout = "{}/{}/f{:0>5}_e{:0>4}.txt".format(trajfolder, mdtag, event.frame+1, event.tid)
        if os.path.isfile(ncout):
            continue
        cppin = """nativecontacts :{} :1-{} writecontacts {}/{}/f{:0>5}_e{:0>4}.txt includesolvent
go
quit
""".format(event.watid+1, prot_size, outfolder, mdtag, event.frame+1, event.tid)
        contacts_cppin = "contacts_{}.cppin".format(next(tempfile._get_candidate_names()))
        with open(contacts_cppin, "w") as t:
            t.write(cppin)
        comm = "cpptraj -p {}/{}/f{:0>5}.pdb -y {}/{}/f{:0>5}.pdb -i {}".format(trajfolder,
                mdtag, event.frame+1, trajfolder, mdtag, event.frame+1, contacts_cppin)
        output = subprocess.run(comm.split(), stdout=subprocess.PIPE, universal_newlines=True)
        if os.path.isfile(contacts_cppin):
            os.remove(contacts_cppin)
    print("Finished ncontacts of md:", mdtag)

def get_native_contacts(datafile, trajfolder, outfolder, prot_size):
    with open(datafile, "rb") as fin:
        data = pickle.load(fin)
    mds = list(data.keys())
    mds.sort()
    os.makedirs(outfolder, exist_ok=True)
    for md in mds:
        os.makedirs(os.path.join(outfolder, md), exist_ok=True)
        _ncontacts(md, data[md], trajfolder, outfolder, prot_size)

def _hbonds(mdtag, events, trajfolder, outfolder):
    for event in events:
        hbout = "{}/{}/f{:0>5}_e{:0>4}.txt".format(outfolder, mdtag, event.frame+1, event.tid)
        svout = "{}/{}/f{:0>5}_e{:0>4}_solv.txt".format(outfolder, mdtag, event.frame+1, event.tid)
        if os.path.isfile(hbout) and os.path.isfile(svout):
            continue
        cppin = """hbond event_H out {}/{}/f{:0>5}_e{:0>4}.txt solventdonor :{} solventacceptor :{}@O dist 3.5 solvout {}/{}/f{:0>5}_e{:0>4}_solv.txt
go
""".format(outfolder, mdtag, event.frame+1, event.tid, event.watid+1, event.watid+1, outfolder, mdtag, event.frame+1, event.tid)
        hbonds_cppin = "hbonds_{}.cppin".format(next(tempfile._get_candidate_names()))
        with open(hbonds_cppin, "w") as t:
            t.write(cppin)
        comm = "cpptraj -p {}/{}/f{:0>5}.pdb -y {}/{}/f{:0>5}.pdb -i {}".format(trajfolder,
                mdtag, event.frame+1, trajfolder, mdtag, event.frame+1, hbonds_cppin)
        output = subprocess.run(comm.split(), stdout=subprocess.PIPE, universal_newlines=True)
        if os.path.isfile(os.path.join(outfolder, "hbonds.log")):
            with open(os.path.join(outfolder, "hbonds.log"), "a") as log:
                log.write(output.stdout)
        else:
            with open(os.path.join(outfolder, "hbonds.log"), "w") as log:
                log.write(output.stdout)
        if os.path.isfile(hbonds_cppin):
            os.remove(hbonds_cppin)
    print("Finished hbonds of md:", mdtag)

def get_hbonds(datafile, trajfolder, outfolder):
    with open(datafile, "rb") as fin:
        data = pickle.load(fin)
    mds = list(data.keys())
    mds.sort()
    os.makedirs(outfolder, exist_ok=True)
    for md in mds:
        os.makedirs(os.path.join(outfolder, md), exist_ok=True)
        _hbonds(md, data[md], trajfolder, outfolder)

if __name__ == "__main__":
    # python3 04_native_contacts.py dhaa/events_data.dat dhaa/frames dhaa 293
    # dhaa:293, epx:319, lipase:534, hepx:316
    # if len(argv) == 4:
    #     get_hbonds(argv[1], argv[2], os.path.abspath(argv[3]))
    # # elif len(argv) == 5:
    # #     get_native_contacts(argv[1], argv[2], os.path.join(argv[3], "ncontacts"), int(argv[4]))
    # else:
    #     print("Usage: python3 04_ncontacts_hbonds.py event_datafile frames_folder output_folder number_of_residues")
    #     print("Usage: python3 04_ncontacts_hbonds.py dhaa/events_data.dat dhaa/frames dhaa #aa (dhaa:293, epx:319, lipase:534, hepx:316)")
    #     quit(1)
    epochs = ["1_4", "1_8"]
    models = ["O", "T3", "T4"]
    for epoch in epochs:
        for model in models:
            database_name = f"../databases/database_{epoch}_{model}.dat"
            frame_folder_name = f"../frames/frames_{epoch}_{model}"
            simulation_folder = "../../../md/simulations"
            hbonds_folder = f"../hbonds/hb_{epoch}_{model}"
            get_hbonds(database_name,frame_folder_name,hbonds_folder)