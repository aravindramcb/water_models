#copy prepared input files from labbit workstation to here via scp 

PREP_DIR="/home/aravind/PhD_local/dean/md_preparation/new/OPC"
declare -a DIR_LIST=("1A" "1.4A" "1.8A" "2.4A" "3A")
FILES=("*.inpcrd" "*.parm7")
mkdir OPC 
mkdir 1A 1.4A 1.8A 2.4A 3A
current_dir=$(pwd)
for dir in ${DIR_LIST[@]}
do
	cd ${dir}
	pwd
	mkdir 1 2 3 4 5
	cd ${current_dir}
done

for ((i=0; i<5; i++))
do
	# echo $i
	for subdir in 1 2 3 4 5
	do
		echo "copying ${DIR_LIST[i]}/$subdir/"
		#mkdir ${DIR_LIST[i]}/$subdir/input_files
		# rm ${DIR_LIST[i]}/$subdir/e*
		scp aravind:$PREP_DIR/${DIR_LIST[i]}/$subdir/*.pdb \
		aravind:$PREP_DIR/${DIR_LIST[i]}/$subdir/*.inpcrd \
		aravind:$PREP_DIR/${DIR_LIST[i]}/$subdir/*.parm7 ${DIR_LIST[i]}/$subdir/input_files/
	done
done



