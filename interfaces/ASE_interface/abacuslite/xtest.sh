#!/bin/bash

set -e # Exit on error

python3 core.py
python3 ./io/generalio.py -v
python3 ./io/legacyio.py -v
python3 ./io/latestio.py -v
python3 ./utils/ksampling.py -v

