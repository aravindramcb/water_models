# Get representative tunnels of P1, to use in figure2
from ctypes import Union

import numpy as np
from transport_tools.libs import tools
import os
import logging


logging.basicConfig(filename='p1_search.log', filemode='w', format='%(levelname)s:%(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)
sh = logging.StreamHandler()
sh.setLevel(logging.DEBUG)
logger.addHandler(sh)

def search_p1_rep_tunnels(scid,group,cutoff:float=0.05):
    mol_system = tools.load_checkpoint(os.path.join("/data/aravindramt/dean/tt/tt_0_9_5/", "_internal", "checkpoints",
                                                    "stage010.dump"))

    p1_bottleneck_mean = {'group1': 1.085346, 'group2': 1.167319,
                          'group3': 1.427511, 'group4': 2.077196,
                          'group5': 2.453579}

    p1_bottleneck_std = {'group1': 0.127381,
                         'group2': 0.145224,
                         'group3': 0.292301,
                         'group4': 0.338774,
                         'group5': 0.337932}
    # Length is the length of P1 tunnel. Here i'm just using SCID 1's avg len and curvature
    p1_avg_len=13.61
    p1_std_len=1.941
    p1_avg_cur=1.231
    p1_std_cur=0.093

    def get_min_max(avg: float, StDev: float):
        """
        This returns min and max values by subtracting avg from std
        :param avg:
        :param StDev:
        :return: (min_value, max_value)
        """
        return avg - StDev, avg + StDev

    min_br,max_br =get_min_max(p1_bottleneck_mean[group],p1_bottleneck_std[group])
    min_len,max_len = get_min_max(p1_avg_len,p1_std_len)
    min_cur,max_cur = get_min_max(p1_avg_cur,p1_std_cur)
    logger.info("Creating Filters")
    filters = tools.define_filters(min_bottleneck_radius=min_br,max_bottleneck_radius=max_br,
                                   min_length=min_len,max_length=max_len,
                                   min_curvature=min_cur,max_curvature=max_cur)
    logger.info("Applying Filters")
    dataset_br = mol_system.get_property_time_evolution_data('bottleneck_radius',
                                                             active_filters=filters,
                                                             sc_id=scid)
    parameter_len = 'length'
    dataset_len = mol_system.get_property_time_evolution_data(parameter_len, active_filters=filters, sc_id=scid)
    parameter_cur = 'curvature'
    dataset_cur = mol_system.get_property_time_evolution_data(parameter_cur, active_filters=filters, sc_id=scid)
    logger.info("Searching frames...")
    def _search_frame():
        number_of_hits = 0
        selected: dict[str, Union[int, list[int]]] = {}
        for md_tag in dataset_br[scid]:
            frame_number = 1
            for frame in zip(dataset_br[scid][md_tag], dataset_len[scid][md_tag], dataset_cur[scid][md_tag]):
                bottleneck = np.round_(frame[0],3)
                length = np.round(frame[1], 3)
                cur = np.round(frame[2], 3)
                bl_len_cur = (bottleneck, length, cur)
                bottleneck_difference = (p1_bottleneck_mean[group]-bottleneck) ** 2
                length_difference = (p1_avg_len - length) ** 2
                curvature_difference = (p1_avg_cur - cur) ** 2
                cutoff_squared = cutoff * cutoff
                if (bottleneck_difference <= cutoff_squared) and (length_difference <= cutoff_squared) \
                        and (curvature_difference <= cutoff_squared):
                    logger.debug(f"{md_tag} -Frame {frame_number} has BR - {bottleneck}, Len - {length}"
                                f" and curvature {cur}")
                    if md_tag not in selected:
                        selected[md_tag] = [frame_number]
                    else:
                        selected[md_tag].append(frame_number)
                    number_of_hits += 1
                frame_number += 1
        return selected,number_of_hits
    selected, number_hits = _search_frame()
    logger.info(f"Number of hits is {number_hits}")
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

def visualize(sc_id,frames,md_tags,group,save_location):
    filters = tools.define_filters()
    out = os.path.join(save_location,f"{group}")
    traj = os.path.join("/data/aravindramt/dean/md/simulations",md_tags,"merged.nc")
    top = os.path.join("/data/aravindramt/dean/md/simulations",md_tags,"structure_HMR.parm7")
    mol_system = tools.load_checkpoint(os.path.join("/data/aravindramt/dean/tt/tt_0_9_5/", "_internal", "checkpoints",
                                                    "stage010.dump"))
    mol_system.show_tunnels_passing_filter(sc_id=sc_id,active_filters=filters,snap_id_list=frames,
                                           out_folder_path=out,md_labels=md_tags)


if __name__ == '__main__':
    save_location = "/home/aravind/PhD_local/dean/figures/main_images/tunnels_representation/figure2_rep_tunnels"
    # tunnel_id = 1
    group = "group5"
    # selections = search_p1_rep_tunnels(tunnel_id,group)

    # md_tags = "1A_opc_1"
    # frames=[38,52,157,441]

    # md_tags = "1.4A_opc_1"
    # frames =[365,408,687,849]
    #
    # md_tags = "1.8A_opc_1"
    # frames = [810,2630,5291,6008,6025]

    # md_tags = "2.4A_opc_4"
    # frames=[7090,7407,8120,10627,15876]

    tunnel_id = 2
    md_tags = '3A_opc_1'
    frames = [681,2184,4675,5010,5966]

    visualize(tunnel_id, frames, md_tags, group, save_location)

