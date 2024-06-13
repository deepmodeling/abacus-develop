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
	# delete time.json (if it exists)
	#--------------------------------------------
	time_json_file="$directory/time.json"
	test -e "$time_json_file" && rm $time_json_file

	#--------------------------------------------
	# delete OUT.autotest (if it exists) 
	#--------------------------------------------
	OUT_directory="$directory/OUT.autotest"
	#test -e "$OUT_directory" && echo $OUT_directory
	test -e "$OUT_directory" && rm -rf $OUT_directory

	#--------------------------------------------
	# delete time.json (if it exists) 
	#--------------------------------------------
	time_json="$directory/time.json"
	#test -e "$time_json" && echo $time_json
	test -e "$time_json" && rm $time_json

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
	sumfile4="$directory/sum_ENV_H2_cube"
	test -e "$sumfile1" && rm -rf $sumfile1
	test -e "$sumfile2" && rm -rf $sumfile2
	test -e "$sumfile3" && rm -rf $sumfile3
	test -e "$sumfile4" && rm -rf $sumfile4

	#--------------------------------------------
	# delete KPT file in kspacing test (if it exists)
	#--------------------------------------------
	if [ ${directory} == "121_PW_kspacing" ]; then
		KPT_file="$directory/KPT"
		test -e "$KPT_file" && rm -rf $KPT_file
	fi

	#--------------------------------------------
	# delete redundant files in RPA test (if it exists)
	#--------------------------------------------
	if [ ${directory} == "282_NO_RPA" ]; then
		redundantfile1="$directory/Cs_data.txt"
		redundantfile2="$directory/Hmpi_0"
		redundantfile3="$directory/Hmpi_1"
		redundantfile4="$directory/Hmpi_2"
		redundantfile5="$directory/Hmpi_3"
		redundantfile6="$directory/dm3_0"
		redundantfile7="$directory/dm3_1"
		redundantfile8="$directory/dm3_2"
		redundantfile9="$directory/dm3_3"
		redundantfile10="$directory/time_0"
		redundantfile11="$directory/time_1"
		redundantfile12="$directory/time_2"
		redundantfile13="$directory/time_3"
		redundantfile14="$directory/KS_eigenvector_0.dat"
		redundantfile15="$directory/band_out"
		redundantfile16="$directory/stru_out"
		redundantfile17="$directory/coulomb_mat_0.txt"
		redundantfile18="$directory/coulomb_mat_1.txt"
		redundantfile19="$directory/coulomb_mat_2.txt"
		redundantfile20="$directory/coulomb_mat_3.txt"
		test -e $redundantfile1  && rm -rf $redundantfile1
		test -e $redundantfile2  && rm -rf $redundantfile2
		test -e $redundantfile3  && rm -rf $redundantfile3
		test -e $redundantfile4  && rm -rf $redundantfile4
		test -e $redundantfile5  && rm -rf $redundantfile5
		test -e $redundantfile6  && rm -rf $redundantfile6
		test -e $redundantfile7  && rm -rf $redundantfile7
		test -e $redundantfile8  && rm -rf $redundantfile8
		test -e $redundantfile9  && rm -rf $redundantfile9
		test -e $redundantfile10 && rm -rf $redundantfile10
		test -e $redundantfile11 && rm -rf $redundantfile11
		test -e $redundantfile12 && rm -rf $redundantfile12
		test -e $redundantfile13 && rm -rf $redundantfile13
		test -e $redundantfile14 && rm -rf $redundantfile14
		test -e $redundantfile15 && rm -rf $redundantfile15
		test -e $redundantfile16 && rm -rf $redundantfile16
		test -e $redundantfile17 && rm -rf $redundantfile17
		test -e $redundantfile18 && rm -rf $redundantfile18
		test -e $redundantfile19 && rm -rf $redundantfile19
		test -e $redundantfile20 && rm -rf $redundantfile20
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
