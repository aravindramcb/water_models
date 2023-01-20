MODEL=$1
ROOT_DIR=/run/user/1000/gvfs/sftp:host=gpu/data/aravindramt/dean/md/${MODEL}
cd $ROOT_DIR
#mkdir minimal_data_test_${MODEL}
for DIR in $(ls -d */)
#for DIR in 1*
do
	cd $DIR
	SUBDIR_ROOT=$(pwd)
	for SUBDIR in $(ls -d */)
	do
		# echo $SUBDIR
#		cp -rnv ${SUBDIR}caver_analyses ../../minimal_data_${MODEL}/${DIR:0:-1}_${SUBDIR}/
		mkdir -p ../../minimal_data_${MODEL}/${DIR:0:-1}_${SUBDIR}aquaduct
		cp -rnv ${SUBDIR}aquaduct/5*.txt ../../minimal_data_${MODEL}/${DIR:0:-1}_${SUBDIR}aquaduct/
		cp -rnv ${SUBDIR}aquaduct/6*.tar.gz ../../minimal_data_${MODEL}/${DIR:0:-1}_${SUBDIR}aquaduct/
	done
	cd $ROOT_DIR
done
