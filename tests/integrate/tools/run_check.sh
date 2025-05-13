#!/bin/bash

CATCH_SCRIPT="../../integrate/tools/catch_properties.sh"
GENERAL_INFO_FILE="../../integrate/general_info"

# check_out: checking the output information
# input: result.out
check_out(){
	outfile=$1
	worddir=`awk '{print $1}' $outfile`
	for word in $worddir; do
        # serach result.out and get information
		cal=`grep "$word" $outfile | awk '{printf "%.'$CA'f\n",$2}'`
        # echo "cal = $cal"
        # search result.ref and get information
		ref=`grep "$word" result.ref | awk '{printf "%.'$CA'f\n",$2}'`
        # echo "ref = $ref"
        # compute the error between 'cal' and 'ref'
		error=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$CA'f\n",x-y}'`
        # echo "error = $error"
        # compare the total time
		if [ $word == "totaltimeref" ]; then		
			echo "$word this-test: $cal ref-1core-test: $ref"
			break
		fi
        # except the total time, other comparison may lead to 'wrong' results
		if [ $(echo "$error == 0"|bc) = 0 ]; then
			echo "----------Wrong!----------"
			echo "word cal ref error"
			echo "$word $cal $ref $error"
			break
		fi
	done
}

test -e $GENERAL_INFO_FILE|| echo "current dir:`pwd`, plese prepare the general_info file."

test -e $GENERAL_INFO_FILE|| exit 0

exec_path=`grep EXEC $GENERAL_INFO_FILE | awk '{printf $2}'`

test -e $exec_path || echo "Error! ABACUS path was wrong!!"

test -e $exec_path || exit 0

CA=`grep CHECKACCURACY $GENERAL_INFO_FILE | awk '{printf $2}'`

NP=`grep NUMBEROFPROCESS $GENERAL_INFO_FILE | awk '{printf $2}'`

path_here=`pwd`
echo "Test in $path_here"

echo "Begin testing the example with $NP cores"

#parallel test
mpirun -np $NP $exec_path > log.txt

test -d OUT.autotest || echo "Some errors occured in ABACUS!"

test -d OUT.autotest || exit 0

#----------------------------------------------------------------------
#if any input parameters for this script, just generate reference file.
#----------------------------------------------------------------------
if test -z $1 
then
$CATCH_SCRIPT result.out
check_out result.out
elif [ $1 == "debug" ] 
then
$CATCH_SCRIPT result.out
check_out result.out
else
$CATCH_SCRIPT result.ref
rm -r OUT.autotest
rm log.txt
fi
