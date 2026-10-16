#!/bin/bash

set -e # Exit on error

python3 ./scf.py -v
python3 ./relax.py -v
python3 ./md.py -v
python3 ./band.py -v
python3 ./magnetic.py -v
if [ "${ASE_ABACUS_SOCKET:-0}" = 1 ]; then
    python3 ./socket_handoff.py -v
fi
