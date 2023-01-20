#!/usr/bin/env bash
# Author <Aravind Thirunavukarasu> <arathi@amu.edu.pl, aravind1233@gmail.com>
# ------- Define Functions -------- #

function tleap_opc {
	     cat > tleap_opc.in << EOF
source leaprc.protein.ff19SB
source leaprc.water.opc
prot = loadpdb $PDB_NAME
check prot
saveamberparm prot TEMP.parm7 TEMP.inpcrd
quit
EOF
      tleap -f tleap_opc.in &> tleap.out
}
function tleap_tip4pew {
      cat > tleap_tip4pew.in << EOF
source leaprc.protein.ff19SB
source leaprc.water.tip4pew
prot = loadpdb $PDB_NAME
check prot
saveamberparm prot TEMP.parm7 TEMP.inpcrd
EOF
      tleap -f tleap_tip4pew.in &>> tleap.out
}

function tleap_tip3p {
    cat > tleap_tip3p.in << EOF
source leaprc.protein.ff19SB
source leaprc.water.tip3p
prot = loadpdb $PDB_NAME
check prot
saveamberparm prot TEMP.parm7 TEMP.inpcrd
quit
EOF
      tleap -f tleap_tip3p.in &>> tleap.out
}

function run_cpptraj {
  cat > cpptraj.in << EOF1
parm TEMP.parm7
parmbox parm TEMP.parm7 truncoct x $DIM y $DIM z $DIM alpha 109.4712190 beta 109.4712190 gamma 109.4712190
parmwrite out TEMP_BOX.parm7
EOF1
    cpptraj -i cpptraj.in &>> cpptraj.log
}

function check_output {
    cat > 05_cpptraj.in << EOF2
parm structure_HMR.parm7
trajin structure.inpcrd
check reportfile check.txt
go
EOF2

    cpptraj 05_cpptraj.in
    grep -v "EPW" check.txt > check_heavy.txt
}

function run_parmed {
    cat > parmed_HMR.in << EOF2
parm TEMP_BOX.parm7
HMassRepartition
parmout structure_HMR.parm7
goX
EOF2
    parmed -O -i parmed_HMR.in
}
# -------- MAIN --------- #
function main {
    cd $ROOT_DIR
    for DIR in $(ls -d */)
    do
      cd $DIR
      SUBDIR_ROOT=$(pwd)
      for SUBDIR in $(ls -d */)
      do
        cd $SUBDIR/temp
        # Copy the prepared pdb to current folder
        PDB_NAME=$(ls e*.pdb)
#        rm $PDB_NAME
#        cp $ORIG_PDB_LOC/$PDB_NAME ./
        # Create new water model based on prevailing water positions
        echo "running tleap"
        tleap_$MODEL_TO_PREPARE
        # Modify topology to truncated octahedral
        DIM=`grep "CRYST1" $PDB_NAME| tr -s [:blank:]| cut -f2 -d" "`
        run_cpptraj
        # modify inpcrd to truncated octahedral
        ChBox -c TEMP.inpcrd -o structure.inpcrd -al 109.4712190 -bt 109.4712190 -gm 109.4712190 -X $DIM -Y $DIM -Z $DIM
    #    rm TEMP_1.inpcrd TEMP_1.parm7 cpptraj.in tleap_tip3p.in
        # Perform HMR
        run_parmed
        # copy final structures to home folder
        cp structure_HMR.parm7 structure.inpcrd ../
        check_output
        echo "Done in $DIR/$SUBDIR"
        cd $SUBDIR_ROOT
      done
      cd $ROOT_DIR
    done
}

# Prepare system for different water models with the extracted PDB from trajectory as input.
MODELS=("tip4pew" "opc" "tip3p")
for ((i=0; i<3; i++))
do
  ORIG_PDB_LOC="/home/aravind/PhD_local/dean/starting_structures/pdb_new/"
  ROOT_DIR="/home/aravind/PhD_local/dean/md_preparation/${MODELS[i]}"
  MODEL_TO_PREPARE=${MODELS[i]}
  main
done