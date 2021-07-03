#!/bin/bash

# ABACUS executable path
abacus=ABACUS
# number of cores
np=2
# threshold with unit: eV
threshold=0.0000001
# check accuracy
ca=8
# regex of case name
case="^[^#].*_PW_.*$"

while getopts a:n:t:c:r: flag
do
    case "${flag}" in
        a) abacus=${OPTARG};;
        n) np=${OPTARG};;
		t) threshold=${OPTARG};;
		c) ca=${OPTARG};;
		r) case=${OPTARG};;
    esac
done

echo "-----AUTO TESTS OF ABACUS ------"
echo "ABACUS path: $abacus";
echo "Number of cores: $np";
echo "Test accuracy: $threshold eV."
echo "Check accuaracy: $ca"
echo "Test cases: $case"
echo "--------------------------------"
echo ""

#----------------------------------------------------------
# define a function named 'check_out'
#----------------------------------------------------------
check_out(){
	#------------------------------------------------------
	# input file $1 is 'result.out' in each test directory
	#------------------------------------------------------
	outfile=$1

	#------------------------------------------------------
	# outfile = result.out
	#------------------------------------------------------
	properties=`awk '{print $1}' $outfile`

	#------------------------------------------------------
	# jd = job description
	#------------------------------------------------------
	if test -e "jd"; then
		jd=`cat jd`
 		echo $jd
	fi

	#------------------------------------------------------
	# check every 'key' word
	#------------------------------------------------------
	for key in $properties; do

		#--------------------------------------------------
		# calculated value
		#--------------------------------------------------
		cal=`grep "$key" $outfile | awk '{printf "%.'$ca'f\n",$2}'`

		#--------------------------------------------------
		# reference value
		#--------------------------------------------------

		ref=`grep "$key" result.ref | awk '{printf "%.'$ca'f\n",$2}'`

		#--------------------------------------------------
		# computed the deviation between the calculated
		# and reference value
		#--------------------------------------------------
		deviation=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$CA'f\n",x-y}'`
		deviation1=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$CA'f\n",y-x}'`

		if [ $key == "totaltimeref" ]; then
			# echo "time=$cal ref=$ref"
			break
		fi


		#--------------------------------------------------
		# If deviation < threshold, then the test passes,
		# otherwise, the test prints out warning
		# Daye Zheng found bug on 2021-06-20,
		# deviation should be positively defined
		#--------------------------------------------------
		if [ $(echo "sqrt($deviation*$deviation) < $threshold"|bc) = 0 ]; then
			echo " *************"
			echo -e "\e[1;31m Incorrect :(  \e[0m"
			echo " *************"
			echo "$key cal=$cal ref=$ref deviation=$deviation"
			break
		else
			#echo "$key cal=$cal ref=$ref deviation=$deviation"
			#echo "[ PASS ] $key"
			echo -e "\e[1;32m [ PASS ] \e[0m $key"
		fi
	done
}

#---------------------------------------------
# the file name that contains all of the tests
#---------------------------------------------

test -e CASES || echo "Plese specify tests"
test -e CASES || exit 0
which $abacus || echo "Error! ABACUS path was wrong!!"
which $abacus || exit 0

testdir=`cat CASES | grep -E $case`

for dir in $testdir; do
	cd $dir
	echo "$dir ($np cores)"
	TIMEFORMAT='Time elapsed: %R seconds'
	#parallel test
	time {
		mpirun -np $np $abacus > log.txt
		test -d OUT.autotest || echo "Some errors happened in ABACUS!"
		test -d OUT.autotest || exit 0

		if test -z $1
		then
			../tools/catch_properties.sh result.out
			check_out result.out
		else
			../tools/catch_properties.sh result.ref
		fi
	}

	echo ""
	cd ../
done
