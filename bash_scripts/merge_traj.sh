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
   		echo "now in $SUBDIR_ROOT/$SUBDIR"
   		
		# Merge the trajectories, image and allign them
    function merge_traj() {
        cat > merge_trajs.sh << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=cpu
#SBATCH --job-name align
#SBATCH --output=%u-slurm_%j.out
#SBATCH --time=2:00:00

module load amber &> /dev/null
module load cuda

if [ "\$SLURM_CLUSTER_NAME" == "labbitgpu" ]
then
    date
    module add amber/20
    export TMPDIR="/tmp/\${USER}/\${SLURM_JOBID}_slurm"
    SLURM_OUTFILE=/tmp/\${USER}-slurm_\${SLURM_JOBID}.out
    mkdir -p \${TMPDIR}
    cd \${TMPDIR}

    echo "------------------------------------------------------"
    echo "Job '\$SLURM_JOBID' was submitted from '\$SLURM_SUBMIT_HOST:\$SLURM_SUBMIT_DIR' to queue '\$SLURM_JOB_PARTITION'."
    echo "Job name is '\$SLURM_JOB_NAME' and is running on the node '\$SLURMD_NODENAME' in temporary folder '\$TMPDIR'"
    echo "------------------------------------------------------"

    FILES_TO_COPY=(structure_HMR.parm7 prod2.nc prod3.nc)
    for COPY_FILE in \${FILES_TO_COPY[*]}
    do
        scp -r \${SLURM_SUBMIT_HOST}:\${SLURM_SUBMIT_DIR}/\${COPY_FILE} \${TMPDIR}/.
    done
function sync_SLURM_output {
scp \${SLURM_OUTFILE} \${SLURM_SUBMIT_HOST}:\${SLURM_SUBMIT_DIR}/slurm_\${SLURM_JOBID}.out
rm \${SLURM_OUTFILE}
}



function sync_file {
FILE_TO_SYNC="\$1"

echo "Syncing file \$FILE_TO_SYNC to \$SLURM_SUBMIT_HOST:\$SLURM_SUBMIT_DIR"
if [ "\$SLURM_CLUSTER_NAME" == "labbitgpu" ]
then
    scp -r \${TMPDIR}/\${FILE_TO_SYNC} \${SLURM_SUBMIT_HOST}:\${SLURM_SUBMIT_DIR}/.
    rm -r \${FILE_TO_SYNC}
fi
}

cat > cpptraj.in << EOF1
parm structure_HMR.parm7
trajin prod/prod2.nc
trajin prod/prod3.nc
autoimage origin familiar
rmsd @N,CA,C mass first out rmsd.dat
atomicfluct out rmsf.dat byres
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

rm -rf prod2.nc prod3.nc cpptraj.in

sync_file "*"

function clean_now {
echo "Cleaning the temporary folder $TMPDIR"
date
sync_SLURM_output
rm -r \${TMPDIR}
}
if [ "\$SLURM_CLUSTER_NAME" == "labbitgpu" ]
then
    trap clean_now INT TERM EXIT SIGINT SIGTERM SIGHUP SIGQUIT SIGABRT
fi

EOF
		chmod u+x merge_trajs.sh
        sbatch merge_trajs.sh
    }
		if [ -f out/prod3.out ]; then
		  merge_traj
		# elif [ ! -f out/prod3.out ]; then
		#   echo "Production run missing"
		#   sbatch 03_production_dean.sh
		fi
		
		cd $SUBDIR_ROOT
	done
	cd $ROOT_DIR
done
