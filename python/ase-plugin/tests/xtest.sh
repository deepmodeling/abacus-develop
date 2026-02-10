#!/bin/bash

set -e # Exit on error

python3 ./scf.py -v
python3 ./relax.py -v
python3 ./md.py -v
python3 ./band.py -v
