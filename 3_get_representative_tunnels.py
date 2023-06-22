# -*- coding: utf-8 -*-

# Get representative tunnel clusters for tt cluster, adopted from custom analysis in TransportTools
__author__ = 'Aravind Selvaram Thirunavukarasu'
__email__ = 'aravind1233@gmail.com'

import logging
import os
from logging import getLogger
from typing import Dict, Any, Union, List

import numpy as np
import pandas as pd

from transport_tools.libs import tools
from transport_tools.libs.config import AnalysisConfig

logging.basicConfig(filename='custom.log', filemode='a', format='%(levelname)s:%(message)s', level=logging.INFO)
logger = getLogger(__name__)
config = AnalysisConfig("/data/aravindramt/dean/tt/config.in")
# load outcomes from tt_engine run
mol_system = tools.load_checkpoint(os.path.join(config.get_parameter("output_path"), "_internal", "checkpoints",
                                                "stage010.dump"))


def print_details(scid, cutoff: float = 0.05) -> dict[str, List[int]]:
    f"""
    Finds the frame where it is +- near the {cutoff} value for Bottleneck and Length for given SC_ID
    :param cutoff: Max deviation of value from actual bottleneck radii and length 
    :param scid: Super cluster ID from transport tools results
    :return: the selected MD_tag and frame numbers
    """
    root_dir = "/data/aravindramt/dean/tt/tt_0_9_5"
    df = pd.read_csv(os.path.join(root_dir + '/statistics/4-filtered_events_statistics.txt'),
                     index_col=False, skiprows=20, skipfooter=1, engine='python', skipinitialspace=True, na_values=0)

    def get_values(dataframe_result, sc_id):
        """
        This returns the avg and stDev of length and bottleneck for requested SC_ID
        """
        _df = dataframe_result.loc[(dataframe_result.SC_ID == sc_id)]
        _Avg_BR = _df.Avg_BR.values[0]
        _StDev_BR = _df.StDev.values[0]
        _Avg_Len = _df.Avg_Len.values[0]
        _StDev_Len = _df["StDev.1"].values[0]
        _Avg_cur = _df.Avg_Cur.values[0]
        _StDev_cur = _df["StDev.2"].values[0]
        return _Avg_BR, _StDev_BR, _Avg_Len, _StDev_Len, _Avg_cur, _StDev_cur

    def get_min_max(avg: float, StDev: float):
        """
        This returns min and max values by subtracting avg from std
        :param avg:
        :param StDev:
        :return: (min_value, max_value)
        """
        return avg - StDev, avg + StDev

    avg_br, std_br, avg_len, std_len, avg_cur, std_cur = get_values(df, sc_id=scid)
    min_len, max_len = get_min_max(avg_len, std_len)
    min_br, max_br = get_min_max(avg_br, std_br)
    min_cur, max_cur = get_min_max(avg_cur, std_cur)

    # What is the length, bottleneck radii of given SC-ID for given snapshots
    # Required - min_length, max_length, min_bottleneck_radius, max_bottleneck_radius, min_num_sims
    filters = tools.define_filters(min_length=min_len, min_bottleneck_radius=min_br,
                                   max_length=max_len, max_bottleneck_radius=max_br
                                   )
    param_br = 'bottleneck_radius'
    dataset_br = mol_system.get_property_time_evolution_data(param_br, active_filters=filters, sc_id=scid)
    parameter_len = 'length'
    dataset_len = mol_system.get_property_time_evolution_data(parameter_len, active_filters=filters, sc_id=scid)
    parameter_cur = 'curvature'
    dataset_cur = mol_system.get_property_time_evolution_data(parameter_cur, active_filters=filters, sc_id=scid)
    logger.info("**************************************************************************************")
    logger.info(f"Average Bottleneck for SCID {scid} is {avg_br}")
    logger.info(f"Average length for SCID {scid} is {avg_len}")
    logger.info(f"Average curvature for {scid} is {avg_cur}")
    logger.info(" ")

    def _search_frame():
        number_of_hits = 0
        selected: dict[str, Union[int, list[int]]] = {}
        for md_tag in dataset_br[scid].keys():
            logger.info(f"Now searching in {md_tag}")
            frame_number = 1
            for frame in zip(dataset_br[scid][md_tag], dataset_len[scid][md_tag], dataset_cur[scid][md_tag]):
                bottleneck = np.round(frame[0], 3)
                length = np.round(frame[1], 3)
                cur = np.round(frame[2], 3)
                bl_len_cur = (bottleneck, length, cur)
                len(bl_len_cur)
                bottleneck_difference = (avg_br - bottleneck) ** 2
                length_difference = (avg_len - length) ** 2
                curvature_difference = (avg_cur - cur) ** 2
                # 0.05 ^ 2 = 0.0025, 0.1 ^ 2 = 0.01 , 0.5 ^ 2 = 0.25
                cutoff_squared = cutoff * cutoff
                if (bottleneck_difference <= cutoff_squared) and (length_difference <= cutoff_squared) \
                        and (curvature_difference <= cutoff_squared):
                    logger.info(f"{md_tag} -Frame {frame_number} has BR - {bottleneck}, Len - {length}"
                                f" and curvature {cur}")
                    if md_tag not in selected:
                        selected[md_tag] = [frame_number]  # See if key exist in selected dict, if not create one.
                    else:
                        selected[md_tag].append(frame_number)
                    number_of_hits += 1
                frame_number += 1
            if number_of_hits == 0:
                logger.info("No nearby values found, moving to next MD_tag")
                continue
            if number_of_hits < 5:
                logger.info("Number of hits less than 5, searching in next MD-Tag")
                continue
            # else:
            #     break
        return selected, number_of_hits

    selected, number_hits = _search_frame()

    while not number_hits >= 5:
        if number_hits != 0 and number_hits >= 5:
            logger.info(f"Number hits is {number_hits}")
            break
        else:
            cutoff = cutoff + 0.05

            logger.info(f"None found with current cutoff of 0.05, extending it by 0.05 and searching with {cutoff}")
            selected, number_hits = _search_frame()
            print(f"NUMBER OF HITS IS {number_hits}")
    return selected


def save_representative_tunnel(scid: int, frame_numbers: list[int], md_labels: list[str], save_dir: str):
    """
    Saves the representative tunnel cluster
    :param save_dir: The name of directory to save the visualization script
    :param scid: Super Cluster ID
    :param frame_numbers: Frame numbers to save
    :param md_labels: MD_tags
    :return: Saves pymol visualization scripts
    """

    out = os.path.join("custom_analysis", save_dir)
    filters = tools.define_filters()
    for md_tag, frames in zip(md_labels, frame_numbers):
        logger.info(f"------Visualization of supercluster {scid} - frames {frames} - {md_tag}---------")
        # noinspection PyTypeChecker
        # modified tools.py to append
        mol_system.show_tunnels_passing_filter(scid, filters, out, md_labels=[md_tag],
                                               snap_id_list=frames, trajectory=False)


def find_frame_and_save_visualization(sc_ids: list[int], save_dir_name: str = 'rep_clusters'):
    sim_frames: dict[int, dict[str, Any]] = {}
    for scid in sc_ids:
        selected = print_details(scid)
        sim_frames[scid] = selected
    for scid in sim_frames:
        # md_tag = list(sim_frames[scid].keys())[0]  # First md_tag
        md_tags = list(sim_frames[scid].keys())
        # frame_number = sim_frames[scid][md_tag][0]  # First value of the md_tag
        frame_numbers = list(sim_frames[scid].values())
        save_representative_tunnel(scid, frame_numbers, md_tags, f"{save_dir_name}/{scid}")



if __name__ == '__main__':
    # Visualize representative superclusters for P1, P2 and P3 tunnels
    P1_scs = [1, 2, 5, 7, 12, 30, 31, 42]
    P2_scs = [3, 4, 6, 11, 16, 25, 27, 41, 44, 43, 50, 58]
    P3_scs = [8, 9, 10, 24]

    # Automated saving the first match
    # Usage - find_frame_and_save_visualization([list of sc ids], save_dir_name)
    find_frame_and_save_visualization(P1_scs, "rep_tunnels/P1")

    # MANUAL SEARCH FOR FRAME NUMBER AND MD_TAG FOR GIVEN SCID:
    # result = print_details(SC_ID, cutoff = new_cutoff)
    result = print_details(11, cutoff=0.5)

    # MANUAL saving of visualization
    # super_cluster_id = 17
    # frame_number = 2454
    # md_tag = "3A_opc_1"
    # save_dir_name = "P2"
    # save_representative_tunnel(super_cluster_id, frame_number, md_tag, save_dir_name)

