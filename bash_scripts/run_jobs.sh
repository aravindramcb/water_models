ROOT_DIR="/data/aravindramt/dean/tip4pew"
cd $ROOT_DIR
for DIR in $(ls -d */)
do
	cd $DIR
	SUBDIR_ROOT=$(pwd)
	for SUBDIR in $(ls -d */)
	do
		# chdir to radius/subdir
		cd $SUBDIR
		# rm -rf temp
   		# echo "now in $SUBDIR_ROOT/$SUBDIR"
   		# copy the scripts
   		# cp ~/dean/scripts/01-equi-prod/* ./ 
   		# cp /home/aravind/PhD_local/dean/md_preparation/scripts/prepare_sys.sh ./

		# 1. Minimization and heating on CPU
		# chdir to radius
		#cp /data/aravindramt/dean/scripts/01-equi-prod/01*.sh ./
#		chmod u+x 01_run_minimization_and_heating_on_CPU.sh
		# sbatch 01_run_minimization_and_heating_on_CPU.sh
		
		# 2. Equllibration on GPU
		 # if [ ! -f "emd0.rst" ]; then
		 # 	  echo "Step 1 is incomplete in $DIR$SUBDIR"
			#   sbatch 01_run_minimization_and_heating_on_CPU.sh
		 #      touch equil_submitted.info
		 # elif [ -f "emd0.rst" ]; then
		 #   echo "Minimization and heating is complete, running equilibration on GPU"
		 #   sbatch 02_run_equllibration_on_GPU_dean.sh
		 # fi
		# check previouis step and run next step 
# 		 if [ ! -f "emd2.rst" ]; then
# 		 	# cp /data/aravindramt/dean/scripts/01-equi-prod/02*.sh ./
# 		 	echo "Step2 emd2.rst is absent ${DIR}${SUBDIR}"

# 		 	sbatch 02_run_equllibration_on_GPU_dean.sh
# 		 elif [[ -f "emd2.rst" ]]; then
# #		  	cp /data/aravindramt/dean/scripts/01-equi-prod/03_production_dean.sh ./
# 		  	echo "01 and 02 steps are complete, running production"
# 		  	sbatch 03_production_dean.sh 
# 		 fi
		# rm -rf caver_analyses *.out *.sh 
		rm starting.inpcrd *.in
		mv structure* temp/
		mkdir min_eq prod
		mv emin* min_eq/
		mv emd* min_eq/
		mv prod* prod/

		
		# Merge the trajectories, image and allign them
    function merge_traj() {
        cat > merge_trajs.sh << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition=interactive
#SBATCH --job-name align
#SBATCH --output=%u-slurm_%j.out
#SBATCH --time=2:00:00

module load amber &> /dev/null

cat > cpptraj.in << EOF1
parm structure_HMR.parm7
trajin prod2.nc
trajin prod3.nc
center @CA mass origin
rmsd @N,CA,C mass first
image origin center familiar
trajout merged.nc netcdf
run
EOF1

cpptraj.MPI -i cpptraj.in

cat > get_pdb.in <<EOF2
parm structure_HMR.parm7
trajin merged.nc 1 1
trajout structure.pdb pdb include_ep
run
EOF2

cpptraj -i get_pdb.in

EOF
		sbatch merge_trajs.sh
    }
		# if [ -f out/prod3.out ]; then
		#   merge_traj
		# elif [ ! -f out/prod3.out ]; then
		#   echo "Production run missing"
		#   sbatch 03_production_dean.sh
		# fi
		# if [[ ! -f aquaduct/6_visualize_results.py ]]; then
		# 	echo "$DIR$SUBDIR"
		# fi
		cd $SUBDIR_ROOT
	done
	cd $ROOT_DIR
done
