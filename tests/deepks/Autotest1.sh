#!/bin/bash

# deepks_test executable path
deepks_test=test_deepks
# regex for test cases
case="^[^#].*_.*$"

while getopts a:r flag
do
    case "${flag}" in
        a) deepks_test=${OPTARG};;
	r) case=${OPTARG};;
    esac
done

echo "-----AUTO TESTS OF ABACUS-DEEPKS 2------"
echo "deepks_test path: $deepks_test";
echo "Test cases: $case"
echo "--------------------------------"
echo ""

testdir=`cat CASES1 | grep -E $case`
failed=0

for dir in $testdir
do
	cd $dir
	echo -e "\e[1;32m [  RUN     ]\e[0m $dir"
	echo -e " [  ------  ] test module_deepks components"
	$deepks_test
	state=`echo $?`
	if [ $state != "0" ]; then
		let failed++
		running_path=`echo "./running.log"`
                cat $running_path
	fi
	cd ..
	echo""
done

if [ $failed -eq 0 ]
then
	exit 0
else
	exit 1
fi


