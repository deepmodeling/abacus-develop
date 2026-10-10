#!/bin/bash

TOOLS_DIR="../../integrate/validation_tools/"

# The note for using the script.
# input parameter: 
	#'none': run and check for this example;
	#'debug': copy scripts to debug, then run and check this example;
	#'other': run this example and generate reference file.

if test -z $1
then
	echo "Run this example and check!"
	$TOOLS_DIR/run_check.sh
elif [ $1 == "debug" ]
then
	echo "Begin debug!"
	cp $TOOLS_DIR/run_check.sh ./
	# catch_properties.sh sources props_*.sh and the helper tools from its
	# own directory, so copy the whole directory, not the entry script alone.
	rm -rf ./validation_tools
	cp -r $TOOLS_DIR ./validation_tools
	./run_check.sh debug
else
	echo "Generate file result.ref ."
	$TOOLS_DIR/run_check.sh $1
fi
