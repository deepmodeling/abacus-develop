#!/bin/bash

#input parameter: 
	#none: run and check for this example;
	#debug: copy the scripts to debug, then run and check this example;
	#other: run this example and generate reference file.

if test -z $1
then
	echo "Run this example and check!"
	./../tools/run_check.sh
elif [ $1 == "debug" ]
then
	echo "Begin debug!"
	cp ../tools/run_check.sh ./
	cp ../tools/catch_properties.sh ./
	./run_check.sh debug
else
	echo "Generate file result.ref ."
	./../tools/run_check.sh $1
fi
