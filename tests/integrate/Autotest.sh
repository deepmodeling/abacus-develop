#!/bin/bash

# ABACUS executable path
abacus=abacus
# number of mpi processes
np=4
# threshold with unit: eV
threshold=0.0000001
# check accuracy
ca=8
# regex of case name
case="^[^#].*_.*$"
# enable AddressSanitizer
sanitize=false

while getopts a:n:t:c:s:r:g flag
do
    case "${flag}" in
        a) abacus=${OPTARG};;
        n) np=${OPTARG};;
		t) threshold=${OPTARG};;
		c) ca=${OPTARG};;
		s) sanitize=${OPTARG};;
		r) case=${OPTARG};;
		g) g=true;; #generate test reference
    esac
done

echo "-----AUTO TESTS OF ABACUS ------"
echo "ABACUS path: $abacus";
echo "Number of cores: $np";
echo "Test accuracy: $threshold eV"
echo "Check accuaracy: $ca"
echo "Test cases: $case"
echo "Generate reference: $g"
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
 		echo " [  ------  ] $jd"
	fi

	#------------------------------------------------------
	# check every 'key' word
	#------------------------------------------------------
	for key in $properties; do

		#--------------------------------------------------
		# calculated value
		#--------------------------------------------------
		cal=`grep "$key" result.out | awk '{printf "%.'$ca'f\n",$2}'`

		#--------------------------------------------------
		# reference value
		#--------------------------------------------------
		ref=`grep "$key" result.ref | awk '{printf "%.'$ca'f\n",$2}'`

		#--------------------------------------------------
		# computed the deviation between the calculated
		# and reference value
		#--------------------------------------------------
		deviation=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$ca'f\n",x-y}'`

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
		if [ ! -n "$deviation" ]; then
            echo -e "\e[1;31m [  FAILED  ]  Fatal Error: key $key not found in output. \e[0m"
			let failed++
			failed_case_list+=$dir
			break
        else
			if [ $(echo "sqrt($deviation*$deviation) < $threshold"|bc) = 0 ]; then
				echo -e "\e[1;33m [  FAILED  ] \e[0m"\
					"$key cal=$cal ref=$ref deviation=$deviation"
				let failed++
				failed_case_list+=$dir
				if [ $(echo "sqrt($deviation*$deviation) < $fatal_threshold"|bc) = 0 ]; then
					let fatal++
					fatal_case_list+=$dir
					echo -e "\e[1;31m [  FATAL   ] \e[0m"\
						"An unacceptable deviation occurs."
					calculation=`grep calculation INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
					running_path=`echo "OUT.autotest/running_$calculation"".log"`
					cat $running_path
				fi
				break
			else
				#echo "$key cal=$cal ref=$ref deviation=$deviation"
				#echo "[ PASS ] $key"
				echo -e "\e[1;32m [      OK  ] \e[0m $key"
			fi
		fi
		let ok++
	done
}

#---------------------------------------------
# the file name that contains all of the tests
#---------------------------------------------

test -e CASES || (echo "Plese specify tests." && exit 1)
which $abacus > /dev/null || (echo "No ABACUS executable was found." && exit 1)

testdir=`cat CASES | grep -E $case`
failed=0
failed_case_list=()
ok=0
fatal=0
fatal_case_list=()
fatal_threshold=1
report=""
repo="$(realpath ..)/"

if [ "$sanitize" == true ]; then
	echo "Testing with Address Sanitizer..."
	mkdir ../html
	echo -e "# Address Sanitizer Diagnostics\n" > ../html/README.md
	report=$(realpath ../html/README.md)
fi

for dir in $testdir; do
	cd $dir
	echo -e "\e[1;32m [  RUN     ]\e[0m $dir"
	TIMEFORMAT=' [  ------  ] Time elapsed: %R seconds'
	#parallel test
	time {
		if [ "$sanitize" == true ]; then
			ASAN_OPTIONS="log_path=asan" mpirun -np $np $abacus > log.txt
			echo -e "## Test case ${dir}\n" >> ${report}
			for diagnostic in asan.*; do
				echo -e "### On process id ${diagnostic}\n" >> ${report}
				echo -e "\`\`\`bash" >> ${report}
				cat ${diagnostic} >> ${report}
				echo -e "\`\`\`\n" >> ${report}
			done
		else
			mpirun -np $np $abacus > log.txt
		fi
		#$abacus > log.txt
		test -d OUT.autotest || (echo "No 'OUT.autotest' dir presented. Some errors may happened in ABACUS." && exit 1)
		if test -z $g
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

if [ "$sanitize" == true ]; then
	if [[ `uname` == "Darwin" ]]; then
		sed -i '' "s,${repo},,g" ${report}
	else
		sed -i "s,${repo},,g" ${report}
	fi
fi

if [ -z $g ]
then
if [ $failed -eq 0 ]
then
	echo -e "\e[1;32m [  PASSED  ] \e[0m $ok test cases passed."
else
	echo -e "\e[1;33m [  FAILED  ] \e[0m $failed test cases out of $[ $failed + $ok ] failed."
	echo $failed_case_list
	if [ $fatal -gt 0 ]
	then
		echo -e "\e[1;31m [  FAILED  ] \e[0m $fatal test cases out of $[ $failed + $ok ] produced fatal error."
		echo $fatal_case_list
		exit 1
	fi
fi
else
echo "Generate test cases complete."
fi
