#!/bin/bash
# Script to get specific frames of PDB file from trajectory
EPOCH_HOME="/mnt/HDD/water_transport/suplementary/4e46/htmd_dihedral/run_htmd/data/"
EPOCH=("e6s1_e5s1p0f1660" "e10s1_e9s3p0f1600" "e1s5_conf1" "e3s3_e2s2p0f1430" "e2s5_e1s2p0f1850" "e9s1_e8s2p0f8840" "e4s3_e3s5p0f3540" "e1s3_conf1" "e8s2_e7s5p0f8080" "e8s4_e7s5p0f9500" "e1s3_conf1" "e9s5_e8s1p0f4660" "e1s3_conf1" "e5s4_e3s5p0f1180" "e5s5_e1s2p0f2860" "e5s4_e3s5p0f1180" "e9s5_e8s1p0f4660" "e5s4_e3s5p0f1180" "e9s5_e8s1p0f4660" "e6s2_e5s1p0f620" "e5s4_e3s5p0f1180" "e5s4_e3s5p0f1180" "e6s2_e5s1p0f620" "e5s4_e3s5p0f1180" "e6s2_e5s1p0f620")
SNAPSHOT=(7090 6486 2510 5456 7896 909 9986 6703 8191 3824 7833 7810 4855 8347 4468 9434 6776 8137 6178 4927 7716 7722 4923 8097 4925 7716 7722 4923 8097 4925)
len=${#EPOCH[@]}
for ((i=0; i<$len; i++))
do
  EPOCH_ID=${EPOCH[i]}
  TRAJECTORY=$EPOCH_HOME${EPOCH[i]}"/trajectory.nc"
  PARM=$EPOCH_HOME${EPOCH[i]}"/structure.parm7"
  PDB_DIR="/home/aravind/PhD_local/dean/starting_structures/pdb_new/"
  cat > get_pdb.in <<EOF
  parm $PARM
  trajin $TRAJECTORY ${SNAPSHOT[i]} ${SNAPSHOT[i]}
  autoimage
  trajout $PDB_DIR${EPOCH_ID:0:4}_${SNAPSHOT[i]}.pdb
  trajout $PDB_DIR${EPOCH_ID:0:4}_${SNAPSHOT[i]}.rst
  go
EOF
  cpptraj -i get_pdb.in

done