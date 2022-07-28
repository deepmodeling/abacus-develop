#!/bin/bash

while getopts a:r: flag
do
    case "${flag}" in
        a) executable=${OPTARG};;
	r) dir=${OPTARG};;
    esac
done

echo "-----AUTO TESTS OF KINETIC OPERATOR------"
echo "Test path: $executable";
echo "Test directory: $dir"
echo "--------------------------------"
echo ""

echo $dir

failed=0

cd $dir
pwd
echo -e "\e[1;32m [  RUN     ]\e[0m $dir"
echo -e " [  ------  ] test kinetic operator in module_hamilt"
$executable
state=`echo $?`
if [ $state != "0" ]; then
	let failed++
fi
cd ..
echo""

if [ $failed -eq 0 ]
then
	exit 0
else
	exit 1
fi


