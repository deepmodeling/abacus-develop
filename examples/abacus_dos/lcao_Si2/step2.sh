#!/bin/bash

cpus=4
materials=Si2_diamond

echo "Setup INPUT file for non-self-consistent calculations of $materials."
cat > INPUT << EOF
INPUT_PARAMETERS
#Parameters (General)
kpoint_file		KPT
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
out_dos			1
EOF

echo "Setup KPT for non-self-consistent calculations."
cat > KPT << EOF
K_POINTS
0
Gamma
12 12 12 0 0 0
EOF

$ABACUS-path/ABACUS.fp_mpi-v1.0.1 $cpus 

echo "nscf done."
