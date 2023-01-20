root_dir="/data/aravindramt/dean/md/$1"
ext_dir="/NAS_data/aravindramt"
for DIR in ("1.4A" "1A" "1.8A" "2.4A" "3A"):
do
  cd $DIR
  SUBDIR_ROOT=$(pwd)
	for SUBDIR in $(ls -d */)
	do
  rsync -rRv temp $ext_dir/opc/
	done
	cd $ROOT_DIR
done