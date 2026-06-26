#!/bin/bash

#bsub -q idle -E ./pre_run.sh -oo %J.log touch 000

if [ $((`date +%s`%20)) -le 10 ]; then
    echo `date +%s` " requeue" >> 001
    exit 1
else
    echo `date +%s` " run" >> 001
    exit 0
fi
