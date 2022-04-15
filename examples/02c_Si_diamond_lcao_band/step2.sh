#!/bin/bash

cpus=4
materials=Si2_diamond

echo "Setup INPUT file for non-self-consistent calculations of $materials."
cat > INPUT << EOF
INPUT_PARAMETERS
#Parameters (General)
kpoint_file		KLINES
pseudo_dir          	./
calculation        	nscf
ntype			1
nbands             	8 
symmetry		0
#Parameters (Methos)
basis_type		lcao
#Parameters (Accuracy)
ecutwfc			50
#Parameters (File)
init_chg		file
out_band		1
EOF

echo "Setup KLINES for non-self-consistent calculations."
cat > KLINES << EOF
K_POINTS
6
Line
0.5 0.0 0.5 20
0.0 0.0 0.0 20
0.5 0.5 0.5 20
0.5 0.25 0.75 20
0.375 0.375 0.75 20
0.0 0.0 0.0 1
EOF

$ABACUS-path/ABACUS.fp_mpi-v1.0.1 $cpus 

echo "nscf done."
