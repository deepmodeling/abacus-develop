#!/bin/bash

for module in "_PW_" "_NO_" "_NP_"; do

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
	test -e "$OUT_directory" && rm -rf $OUT_directory

	#--------------------------------------------
    # delete projected_DM.dat (if it exists) 
    #--------------------------------------------
	projected_DM="$directory/projected_DM.dat"
	test -e "$projected_DM" && rm -rf $projected_DM

done

done

test -e "check_file" && rm -rf check_file
