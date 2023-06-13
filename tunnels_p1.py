import os
import pytraj as pt

def save_frames(md_tags, frames, group, save_location):
    traj = os.path.join("/data/aravindramt/dean/md/simulations",md_tags,"merged.nc")
    top = os.path.join("/data/aravindramt/dean/md/simulations",md_tags,"structure_HMR.parm7")
    sim_loaded=pt.load(traj,top,frames)
    sim_loaded=pt.strip(sim_loaded,"!@CA,C,N")
    out = os.path.join(save_location, f"group{group}/frames.pdb")
    pt.write_traj(out,sim_loaded,format='PDB',overwrite=True,options='model')
    print("done")
if __name__ == '__main__':
    save_location = "/home/aravind/PhD_local/dean/figures/main_images/tunnels_representation/figure2_rep_tunnels"

    # md_tags = "1A_opc_1"
    # frames=[38,52,157,441]
    #
    # md_tags = "1.4A_opc_1"
    # frames =[365,408,687,849]
    #
    # md_tags = "1.8A_opc_1"
    # frames = [810,2630,5291,6008,6025]
    #
    # md_tags = "2.4A_opc_4"
    # frames = [7090, 7407, 8120, 10627, 15876]

    md_tags = '3A_opc_1'
    frames = [681,2184,4675,5010,5966]

    group=5
    save_frames(md_tags, frames, group, save_location)

