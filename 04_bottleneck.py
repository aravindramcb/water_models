# python program to process bottleneck residues
# -*- coding: utf-8 -*-

# Author: Aravind Thirunavukarasu
# Email : arathi@amu.edu.pl, aravind1233@gmail.com

import os
from dataclasses import dataclass, asdict


@dataclass
class Data:
    SC_ID: str
    Total_No_Frames: int
    Bottleneck_residues: dict

def get_bottleneck(results_location: str):
    if os.getcwd() != results_location:
        os.chdir(results_location)
    folders = [d for d in os.listdir() if os.path.isdir(d)]
    bottleneck_results = "2-filtered_tunnels_statistics_bottleneck_residues.txt"

    for folder in folders:
        with open(os.path.join(folder, bottleneck_results), "r") as file:
            lines = file.readlines()
            data_string = ''.join(lines[18:])
            data = Data()
            data_objects = []

            for data_string in data_list:
                sc_id, total_no_frames, *residues = data_string.split(', ')

                bottleneck_residues = {}
                for residue in residues:
                    residue_id, frequency = residue.split(':')
                    bottleneck_residues[int(residue_id)] = float(frequency)

                data_objects.append(Data(int(sc_id), int(total_no_frames), bottleneck_residues))

            for data_object in data_objects:
                print(asdict(data_object))


if __name__ == '__main__':
    comparative_analysis_dir = "/data/aravindramt/dean/tt/tt_0_9_5_bottleneck/statistics/comparative_analysis/"
    get_bottleneck(comparative_analysis_dir)

