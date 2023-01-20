ROOT_DIR="/data/aravindramt/dean/OPC/"
RESTART_ID=1
for DIR in $(ls -d */)
do
	cd $DIR
	SUBDIR_ROOT=$(pwd)
	for SUBDIR in $(ls -d */)
	do
		# cd $SUBDIR
   		
   		rsync -Rv $SUBDIR/merged.nc eagle:/tmp/lustre_shared/aravindram/dean/$DIR
	done
	# cd $ROOT_DIR
done
