#!/bin/bash


check_out(){
	outfile=$1
	worddir=`awk '{print $1}' $outfile`
	for word in $worddir; do
		cal=`grep "$word" $outfile | awk '{printf "%.'$CA'f\n",$2}'`
		ref=`grep "$word" result.ref | awk '{printf "%.'$CA'f\n",$2}'`
		error=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$CA'f\n",x-y}'`
		if [ $word == "totaltimeref" ]; then		
			echo "$word this-test: $cal ref-1core-test: $ref"
			break
		fi
		if [ $(echo "$error == 0"|bc) = 0 ]; then
			echo "----------Wrong!----------"
			echo "word cal ref error"
			echo "$word $cal $ref $error"
			break
		fi
	done
}

test -e ../general_info||echo "plese write the file general_info"
test -e ../general_info||exit 0
exec_path=`grep EXEC ../general_info|awk '{printf $2}'`
test -e $exec_path || echo "Error! ABACUS path was wrong!!"
test -e $exec_path || exit 0
CA=`grep CHECKACCURACY ../general_info | awk '{printf $2}'`
NP=`grep NUMBEROFPROCESS ../general_info | awk '{printf $2}'`

	path_here=`pwd`
	echo "Test in $path_here"

	echo "Begin test with $NP cores"
	#parallel test
	mpirun -np $NP $exec_path > log.txt
	test -d OUT.autotest || echo "Some errors happened in ABACUS!"
	test -d OUT.autotest || exit 0

	#if any input parameters for this script, just generate reference file.
	if test -z $1 
	then
		../tools/catch_properties.sh result.out
		check_out result.out
	elif [ $1 == "debug" ] 
	then
		./catch_properties.sh result.out
		check_out result.out
	else
		../tools/catch_properties.sh result.ref
		rm -r OUT.autotest
		rm log.txt
	fi
