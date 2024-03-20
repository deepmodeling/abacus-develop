#!/bin/bash
for i in STRU_*
do
stru=$(echo $i|cut -d"_" -f2)
echo $stru
echo "stru=$stru"
mkdir SCF-$stru
cd SCF-$stru
pwd
cat > INPUT <<EOF
INPUT_PARAMETERS
#Parameters     (General)
suffix          DIA-50-$stru
calculation     scf
esolver_type    ksdft
pseudo_dir      ../
nbands          45
symmetry        1
cal_force       1
cal_stress      1
kpar            2

#Parameters (Accuracy)
ecutwfc         50
scf_thr         1e-12
scf_nmax        100
basis_type      pw
ks_solver       cg
smearing_method gauss
smearing_sigma  0.01
mixing_type     broyden
mixing_beta     0.7

stru_file       STRU_$stru
EOF
cat > KPT <<EOF
K_POINTS
0
Gamma
2 2 2 0 0 0
EOF
cp ../STRU_$stru .
mpirun -n 96 ABACUS.mpi
# sbatch ../sub.sh
cd ../
done
