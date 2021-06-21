#!/bin/bash

# threshold with unit: eV
threshold="0.0000001"

echo "-----AUTO TESTS OF ABACUS ------"
echo "test accuracy is $threshold eV."
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
		cal=`grep "$key" $outfile | awk '{printf "%.'$CA'f\n",$2}'`

		#--------------------------------------------------
		# reference value
		#--------------------------------------------------

		ref=`grep "$key" result.ref | awk '{printf "%.'$CA'f\n",$2}'`
		#--------------------------------------------------
		# computed the deviation between the calculated 
		# and reference value
		#--------------------------------------------------
		deviation=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$CA'f\n",x-y}'`
		deviation1=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$CA'f\n",y-x}'`

		if [ $key == "totaltimeref" ]; then		
#			echo "time=$cal ref=$ref"
			break
		fi

		
		#--------------------------------------------------
		# If deviation < threshold, then the test passes,
		# otherwise, the test prints out warning
		#--------------------------------------------------
		if [ $(echo "$deviation < $threshold && $deviation1 < $threshold"|bc) = 1 ]; then
			#echo "$key cal=$cal ref=$ref deviation=$deviation"
			#echo "[ PASS ] $key"
			echo -e "\e[1;32m [ PASS ] \e[0m $key"
		else 
			echo " *************"
			echo -e "\e[1;31m Incorrect :(  \e[0m"
			echo " *************"
			echo "$key cal=$cal ref=$ref deviation=$deviation"
			break
		fi
	done
}

#---------------------------------------------
# the file name that contains all of the tests
#---------------------------------------------

test -e general_info||echo "plese write the file list_of_tests"
test -e general_info||exit 0
exec_path=`grep EXEC general_info|awk '{printf $2}'`
test -e $exec_path || echo "Error! ABACUS path was wrong!!"
test -e $exec_path || exit 0
CA=`grep CHECKACCURACY general_info | awk '{printf $2}'`
NP=`grep NUMBEROFPROCESS general_info | awk '{printf $2}'`

grep -w TESTDIR general_info | awk '{print $2}' > testdir.txt
list=`grep list_of_tests general_info|sed -e 's/[^ ]* //'`
if [ -z "$list" ]
then
#echo "do no thing"
testdir=`cat testdir.txt`
else
#echo $list
value_line=(` echo $list | head -n1 `)

colume=`echo ${#value_line[@]}`
for (( col=0 ; col<$colume ; col++ ));do
    value=`echo ${value_line[$col]}`
    grep $value testdir.txt > testdir.dat
    mv testdir.dat testdir.txt
done
testdir=`cat testdir.txt`
fi
rm testdir.txt

for dir in $testdir; do
cd $dir
	echo "$dir ($NP cores)"
	#parallel test
	mpirun -np $NP $exec_path > log.txt
	test -d OUT.autotest || echo "Some errors happened in ABACUS!"
	test -d OUT.autotest || exit 0

	if test -z $1 
	then
		../tools/catch_properties.sh result.out
		check_out result.out
	else
		../tools/catch_properties.sh result.ref
	fi

	echo ""
cd ../
done	
