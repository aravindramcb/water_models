ROOT_DIR="/data/aravindramt/dean/tip4pew/"
RESTART_ID=1
cd $ROOT_DIR
for DIR in $(ls -d */)
do
	cd $DIR
	SUBDIR_ROOT=$(pwd)
	for SUBDIR in $(ls -d */)
	do
		cd $SUBDIR
		rm -rf aquaduct
		mkdir aquaduct
		if [ ! -f merged.nc  ]; then
			echo "merged.nc not found"
			exit
		fi
		cd aquaduct
		if [[ ! -f 6_visualize_results.py ]]; then
			echo "Now in $(pwd)"
			
			cat > config.txt <<EOF2
[global]
top = ../structure.pdb
trj = ../merged.nc 

[traceable_residues]
scope = backbone
dump = 1_traceable_residues_data.dump
scope_convexhull = True
object = (resname WAT) and (sphzone 6.0 ((resnum 103 and name CG) or (resnum 104 and name CD2) or (resnum 243 and name N) or (resnum 165 and name N)))

[raw_paths]
clear_in_object_info = False
dump = 2_raw_paths_data.dump

[separate_paths]
dump = 3_separate_paths_data.dump
sort_by_id = True
apply_smoothing = False
apply_soft_smoothing = False
discard_short_paths = 1

[inlets_clusterization]
recluster_outliers = True
dump = 4_inlets_clustering_data.dump
detect_outliers = Auto
singletons_outliers = 2
create_master_paths = True
max_level = 2

[analysis]
save = 5_analysis_results.txt
dump_config = True
scope_chull = backbone
object_chull = resname WAT and (sphzone 6.0 ((resnum 103 and name CG) or (resnum 104 and name CD2) or (resnum 243 and name N) or (resnum 165 and name N)))

[visualize]
all_paths_raw = True
all_paths_smooth = False
all_paths_split = True
all_paths_raw_io = True
all_paths_smooth_io = False
paths_raw = True
paths_smooth = False
paths_states = True
paths_raw_io = True
paths_smooth_io = False
ctypes_raw = True
ctypes_smooth = False
inlets_clusters = True
show_molecule = protein
show_molecule_frames = 0
show_scope_chull =  backbone	
show_object_chull = resname WAT and (sphzone 6.0 ((resnum 103 and name CG) or (resnum 104 and name CD2) or (resnum 243 and name N) or (resnum 165 and name N)))
save = 6_visualize_results.py

[clustering]
method = meanshift
cluster_all = True
bandwidth = Auto
recursive_clusterization = clusterization
recursive_threshold = >0.9

[reclusterization]
method = meanshift
cluster_all = False
bandwidth = Auto
EOF2
			cat > run_aquaduct.sh <<EOF1
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=interactive
#SBATCH --job-name=aq_$DIR$SUBDIR
#SBATCH --mem=8gb
#SBATCH --output=%j.out
			
valve.py -t 4 -c config.txt &> aquaduct_runtime.log
EOF1
		fi
		echo "the directory contains $(dir)"
		chmod u+x run_aquaduct.sh
		sbatch run_aquaduct.sh
		# 	valve.py -t 4 -c config.txt &> aquaduct_runtime.log
		# fi
		cd $SUBDIR_ROOT
	done
	cd $ROOT_DIR
done
