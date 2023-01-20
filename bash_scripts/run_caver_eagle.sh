#!/bin/bash 
ROOT_DIR=$(pwd)
for DIR in $(ls -d */)
do
	cd $DIR
	SUBDIR_ROOT=$(pwd)
	for SUBDIR in $(ls -d */)
	do
		cd $SUBDIR
   	cp ~/grant_312/project_data/aravind/scripts/*.sh ./
   	# 01
   	bash 01-prepare_caver.sh
   	# 02
   	# 03
   	# 05
		cd $SUBDIR_ROOT
	done
	cd $ROOT_DIR
done
