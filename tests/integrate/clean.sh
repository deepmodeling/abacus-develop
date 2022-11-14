#!/bin/bash

for module in "_PW_" "_NO_" "_NP_" "_LJ_" "_OF_" "tools"; do

echo $module

for directory in `ls | grep $module`; do

	#--------------------------------------------
	# delete log.txt (if it exists) 
	#--------------------------------------------
	log_file="$directory/log.txt"
	#test -e "$log_file" && echo $log_file
	test -e "$log_file" && rm $log_file

	#--------------------------------------------
	# delete result.out (if it exists) 
	#--------------------------------------------
	result_file="$directory/result.out"
	#test -e "$result_file" && echo $result_file
	test -e "$result_file" && rm $result_file

	#--------------------------------------------
	# delete OUT.autotest (if it exists) 
	#--------------------------------------------
	OUT_directory="$directory/OUT.autotest"
	#test -e "$OUT_directory" && echo $OUT_directory
	test -e "$OUT_directory" && rm -rf $OUT_directory

	#--------------------------------------------
	# delete descriptor.dat (if it exists) 
	#--------------------------------------------
	descriptor="$directory/descriptor.dat"
	test -e "$descriptor" && rm -rf $descriptor

	#--------------------------------------------
	# delete H_V_delta.dat (if it exists) 
	#--------------------------------------------
	H_V_delta="$directory/H_V_delta.dat"
	test -e "$H_V_delta" && rm -rf $H_V_delta

	#--------------------------------------------
	# delete projected_DM.dat (if it exists) 
	#--------------------------------------------
	projected_DM="$directory/projected_DM.dat"
	test -e "$projected_DM" && rm -rf $projected_DM

	#--------------------------------------------
	# delete Onsager.txt (if it exists)
	#--------------------------------------------
	onsager_file="$directory/Onsager.txt"
	test -e "$onsager_file" && rm -rf $onsager_file

	#--------------------------------------------
	# delete kpoints (if it exists)
	#--------------------------------------------
	kpoints_file="$directory/kpoints"
	test -e "$kpoints_file" && rm -rf $kpoints_file

	#--------------------------------------------
	# delete exec files in tools directory (if it exists)
	#--------------------------------------------
	sumfile1="$directory/sum_BAND_CHG_H2"
	sumfile2="$directory/sum_BAND_CHG_H2_cube"
	sumfile3="$directory/sum_ENV_H2"
	test -e "$sumfile1" && rm -rf $sumfile1
	test -e "$sumfile2" && rm -rf $sumfile2
	test -e "$sumfile3" && rm -rf $sumfile3

	#--------------------------------------------
	# delete KPT file in kspacing test (if it exists)
	#--------------------------------------------
	if [ ${directory} == "121_PW_kspacing" ]; then
		KPT_file="$directory/KPT"
		test -e "$KPT_file" && rm -rf $KPT_file
	fi

	#--------------------------------------------
	# delete *.npy (if it exists) 
	#--------------------------------------------
	num=$(find -name '*.npy' | wc -l)
	if [ $num != "0" ]; then
		rm -rf $directory/*.npy
	fi

	#--------------------------------------------
	# delete test_exx (if it exists) 
	#--------------------------------------------
	test_exx="$directory/test_exx"
	test -e "$test_exx" && rm -rf $test_exx

done

done
