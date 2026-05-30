#!/bin/bash

for ii in *
do
	if [ -d $ii ];then
		cd ${ii}
		echo "RUN: ${ii}"
		bash run.sh
		cd ..
	fi
done
