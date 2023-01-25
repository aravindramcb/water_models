#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2021 Carlos Eduardo Sequeiros Borja <carseq@amu.edu.pl>
#  

import os
import pickle
import argparse
import numpy as np
from sys import argv
import configparser as cfp
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

def _analyze_tevents(infiles, dist_threshold, frac_threshold):
    events = {}
    for f in infiles:
        _id = int(os.path.basename(f).split("_")[0])
        if _id not in events:
            event = TEvent(tid=_id)
            events[_id] = event
        else:
            event = events[_id]
        with open(f, "r") as inf:
            line = inf.readline()
            if "entry" in line:
                event.event_type = "entry"
            else:
                event.event_type = "release"
            _watid = int(line.strip().split("'")[6][4:])
            _sc = int(line.strip().split()[-1])
            line = inf.readline()
            _fraction = float(line.strip().split()[-1])
            if (_fraction >= frac_threshold) and (_fraction > event.fraction):
                [inf.readline() for j in range(4)]
                for line in inf:
                    chunks = [float(c) for c in line.strip().split(",")]
                    if chunks[4] <= dist_threshold and chunks[-2] < event.radius:
                        event.watid = _watid
                        event.supercluster = _sc
                        event.fraction = _fraction
                        event.frame = int(chunks[0])
                        event.radius = chunks[-2]
                        event.dist2lig = chunks[4]
                        event.caver_cluster = int(chunks[5])
    return events

def parse_events(infolder, dist_threshold=0.0, frac_threshold=0.0):
    entries = [os.path.join(infolder, f) for f in os.listdir(infolder) if "entry" in f]
    releases = [os.path.join(infolder, f) for f in os.listdir(infolder) if "release" in f]
    entries.sort()
    releases.sort()
    entries = _analyze_tevents(entries, dist_threshold, frac_threshold)
    releases = _analyze_tevents(releases, dist_threshold, frac_threshold)
    return entries, releases


_sim_list_1_4_opc = ['1.4A_opc_1', '1.4A_opc_2', '1.4A_opc_3', '1.4A_opc_4', '1.4A_opc_5']
_sim_list_1_4_tip3p =[ '1.4A_tip3p_1', '1.4A_tip3p_2','1.4A_tip3p_3', '1.4A_tip3p_4', '1.4A_tip3p_5']
_sim_list_1_4_tip4pew=['1.4A_tip4pew_1', '1.4A_tip4pew_2', '1.4A_tip4pew_3','1.4A_tip4pew_4', '1.4A_tip4pew_5']

_sim_list_1_opc=['1A_opc_1', '1A_opc_2', '1A_opc_3','1A_opc_4', '1A_opc_5']
_sim_list_1_tip3p=[ '1A_tip3p_1', '1A_tip3p_2', '1A_tip3p_3', '1A_tip3p_4', '1A_tip3p_5']
_sim_list_1_tip4pew= [ '1A_tip4pew_1','1A_tip4pew_2', '1A_tip4pew_3', '1A_tip4pew_4', '1A_tip4pew_5']

_sim_list1_8=['1.8A_opc_1', '1.8A_opc_2', '1.8A_opc_3', '1.8A_opc_4', '1.8A_opc_5']
_sim_list1_8_tip3p=['1.8A_tip3p_1', '1.8A_tip3p_2', '1.8A_tip3p_3', '1.8A_tip3p_4', '1.8A_tip3p_5']
_sim_list1_8_tip4pew=[ '1.8A_tip4pew_1','1.8A_tip4pew_2', '1.8A_tip4pew_3', '1.8A_tip4pew_4', '1.8A_tip4pew_5']

others = [  '2.4A_opc_1', '2.4A_opc_2', '2.4A_opc_3',
            '2.4A_opc_4', '2.4A_opc_5', '2.4A_tip3p_1', '2.4A_tip3p_2', '2.4A_tip3p_3', '2.4A_tip3p_4', '2.4A_tip3p_5',
            '2.4A_tip4pew_1', '2.4A_tip4pew_2', '2.4A_tip4pew_3', '2.4A_tip4pew_4', '2.4A_tip4pew_5', '3A_opc_1',
            '3A_opc_2', '3A_opc_3', '3A_opc_4', '3A_opc_5', '3A_tip3p_1', '3A_tip3p_2', '3A_tip3p_3', '3A_tip3p_4',
            '3A_tip3p_5', '3A_tip4pew_1', '3A_tip4pew_2', '3A_tip4pew_3', '3A_tip4pew_4', '3A_tip4pew_5']

def make_event_database():
    # ttconfig = cfp.ConfigParser()
    # ttconfig.read(args.tt_config_ini)
    # tt_path = os.path.dirname(args.tt_config_ini)
    tt_path = '/data/aravindramt/dean/tt/'
    em_path = os.path.join(tt_path, 'tt_0_9_5', "data", "exact_matching_analysis")
    # mds = [f for f in os.listdir(em_path)]
    mds=_sim_list_1_opc
    mds.sort()
    frames_bymd = defaultdict(list)
    for md in mds:
        entries, releases = parse_events(os.path.join(em_path, md), 0.0, 0.7)
        for entry in entries.values():
            if entry.frame > -1:
                frames_bymd[md].append(entry)
        for release in releases.values():
            if release.frame > -1:
                frames_bymd[md].append(release)
    with open('database_1_O.dat', "wb") as fout:
        pickle.dump(frames_bymd, fout)
    return frames_bymd

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Build database of transport events")
    # parser.add_argument("-c", "--ttconfig", action="store", required=True, dest="tt_config_ini",
    #     help="The ini file used in the TransportTools run.")
    # parser.add_argument("-o", "--out", action="store", required=True, dest="output",
    #     help="Filename for the output database.")
    # args = parser.parse_args()
    make_event_database()
