#!/bin/bash

for i in P*;do
    rm -rf $i/result.log $i/time.log $i/OUT.* $i/result.out
done

rm -rf *cpu  *kpar *bxyz sum.dat* *.old log
