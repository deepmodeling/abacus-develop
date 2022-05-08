#!/bin/bash

cpus=4
materials=Si2_diamond

echo "Setup INPUT file for self-consistent calculations of $materials."
cat > INPUT << EOF
INPUT_PARAMETERS
#Parameters	(General)
pseudo_dir          		./
calculation			relax
ntype		  	    	1	
nbands		  		8
#Parameters (Methos)
basis_type			lcao
symmetry			0	
#Parameters (Accuracy)
ecutwfc				50
scf_thr				1.0e-7	// about iteration

relax_nmax				50
force_thr_ev			1.0e-3
#Parameters (File)
out_chg			1
EOF

cat > KPT << EOF
K_POINTS
0
Gamma
4 4 4 0 0 0
EOF

$ABACUS-path/ABACUS.fp_mpi-v1.0.1 $cpus 
echo "scf done."

rm KPT
rm INPUT

#../RunTest 4
