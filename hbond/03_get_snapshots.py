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

def _write_frames(mdtag, events, trajfolder, outfolder):
    frames = set()
    for e in events:
        # We add 1 to the frame described in the exact matching analysis files
        # to correspond to the actual frame in the MD trajectory
        frames.add(e.frame + 1)
    frames = list(frames)
    frames.sort()
    to_run = False
    limit = len(frames)//1000
    i = 0
    while i < limit:
        cppin = """parm {}/structure_HMR.parm7
trajin {}/merged.nc
""".format(os.path.join(trajfolder, mdtag), os.path.join(trajfolder, mdtag))
        for f in frames[i*1000:(i+1)*1000]:
            pdbout = "{}/{}/f{:0>5}.pdb".format(outfolder, mdtag, f)
            if os.path.isfile(pdbout):
                continue
            cppin += "trajout {} onlyframes {}\n".format(pdbout, f)
            to_run = True
        cppin += "go\nquit\n"
        ef_cppin = "extract_frames_{}_{}.cppin".format(mdtag, i)
        with open(ef_cppin, "w") as t:
            t.write(cppin)
        if to_run:
            output = subprocess.run(["cpptraj", "-i", ef_cppin],
                                    stdout=subprocess.PIPE, universal_newlines=True)
            
            if os.path.isfile("extract_frames.log"):
                with open(os.path.join(outfolder, "extract_frames.log"), "a") as log:
                    log.write(output.stdout)
            else:
                with open(os.path.join(outfolder, "extract_frames.log"), "w") as log:
                    log.write(output.stdout)
        i += 1
    
    cppin = """parm {}/structure_HMR.parm7
trajin {}/merged.nc
""".format(os.path.join(trajfolder, mdtag), os.path.join(trajfolder, mdtag))
    for f in frames[i*1000:]:
        pdbout = "{}/{}/f{:0>5}.pdb".format(outfolder, mdtag, f)
        if os.path.isfile(pdbout):
            continue
        cppin += "trajout {} onlyframes {}\n".format(pdbout, f)
        to_run = True
    cppin += "go\nquit\n"
    ef_cppin = "extract_frames_{}_{}.cppin".format(mdtag, i)
    with open(ef_cppin, "w") as t:
        t.write(cppin)
    if to_run:
        output = subprocess.run(["cpptraj", "-i", ef_cppin],
                                stdout=subprocess.PIPE, universal_newlines=True)
        
        if os.path.isfile("extract_frames.log"):
            with open(os.path.join(outfolder, "extract_frames.log"), "a") as log:
                log.write(output.stdout)
        else:
            with open(os.path.join(outfolder, "extract_frames.log"), "w") as log:
                log.write(output.stdout)
    print("Extracted frames for md:", mdtag)

def get_pdb_frames(datafile, trajfolder, outfolder):
    with open(datafile, "rb") as fin:
        data = pickle.load(fin)
    mds = list(data.keys())
    mds.sort()
    os.makedirs(outfolder, exist_ok=True)
    for md in mds:
        os.makedirs(os.path.join(outfolder, md), exist_ok=True)
        _write_frames(md, data[md], trajfolder, outfolder)

if __name__ == "__main__":
    # python3 get_snapshots.py dhaa/events_data.dat dhaa/data dhaa/frames
    # if len(argv) == 4:
    #     get_pdb_frames(argv[1], argv[2], argv[3])
    # else:
    #     print("Usage: python3 03_get_snapshots.py event_datafile trajectory_folders output_folder")
    #     print("Usage: python3 03_get_snapshots.py dhaa/events_data.dat dhaa/data dhaa/frames")
    #     quit(1)
    # epochs = ["1_4", "1_8"]
    epochs =["1"]
    models = ["O", "T3", "T4"]
    for epoch in epochs:
        for model in models:
            database_name = f"../databases/database_{epoch}_{model}.dat"
            frame_folder_name = f"../frames/frames_{epoch}_{model}"
            simulation_folder = "../../../md/simulations"
            get_pdb_frames(database_name,simulation_folder,frame_folder_name)