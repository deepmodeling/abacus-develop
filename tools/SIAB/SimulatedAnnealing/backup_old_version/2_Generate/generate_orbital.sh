#!/bin/bash
# author mohan 
# edit by Pengfei Li 2013-6-4
# edit by Wenshuai Zhang 2016-11-30
#-----------------------------------------------------------------
#
# (0.0) fixed  parameters 
#
#-----------------------------------------------------------------
mass=1                                         #not used yet.
lat0=20                                        #in a.u
cpu_num=8         
export OMP_NUM_THREADS=1
#hostfpath="/home/nic/wszhang/eclipse_project/mesia_dft/error_estimates_for_DFT/cif2cellroot/0abacus_lcao/myhosts"
#hostfpath="/home/nic/wszhang/eclipse_project/mesia_dft/delta_test/delta_dft/cifs2deltaDFT/myhosts"
#-----------------------------------------------------------------
#
# (0.1) input parameters
#
#-----------------------------------------------------------------
# (0.1.1) get exe
EXE_pw=`grep "EXE_pw" ORBITAL_INPUT | awk -F "EXE_pw" '{print $0}' | awk '{print $2}'`

# (0.1.2)get SIA
EXE_orbital=`grep "EXE_orbital" ORBITAL_INPUT | awk -F "EXE_orbital" '{print $0}' | awk '{print $2}'`

# (0.1.3)get the targets element and id
targets=`grep "targets" ORBITAL_INPUT | awk -F "targets" '{print $0}' | awk '{print $2}'`

for name in $targets; do

element=`echo "$name" | awk -F "_" '{print $2}'`
id=`echo "$name" | awk -F "_" '{print $1}'`

# (0.1.4)get the nbands
nbands=`grep "nbands" ORBITAL_INPUT | awk -F "nbands" '{print $0}' | awk '{print $2}'`

# (0.1.5)get the ref_bands
ref_bands=`grep "ref_bands" ORBITAL_INPUT  | awk -F "$ref_bands" '{print $0}' | awk '{print $2}'`

# (0.1.6)get the pseudo_dir
pseudo_dir=`grep "Pseudo_dir" ORBITAL_INPUT | awk -F "Pseudo_dir" '{print $0}' | awk '{print $2}'`

# (0.1.7)get maxL S P D
maxL=`grep "maxL" ORBITAL_INPUT | awk -F "maxL" '{print $0}' | awk '{print $2}'`

#if ( test $maxL = 0 )   // mohan's scheme
#then 
#S_MIN=1
#P_MIN=" "
#D_MIN=" "
#elif ( test $maxL = 1 )
#then
#S_MIN=1
#P_MIN=1
#D_MIN=" "
#else
#S_MIN=1
#P_MIN=1
#D_MIN=1
#fi

# (0.1.8)get the level
Level=`grep "Level" ORBITAL_INPUT | awk -F "level" '{print $0}' | awk '{print $2}'`

# (0.1.9)get every level`s lmax s p d f g
L[0]=`grep "level1" ORBITAL_INPUT | awk -F "level1" '{print $2}'`
L[1]=`grep "level2" ORBITAL_INPUT | awk -F "level2" '{print $2}'`
L[2]=`grep "level3" ORBITAL_INPUT | awk -F "level3" '{print $2}'`
L[3]=`grep "level4" ORBITAL_INPUT | awk -F "level4" '{print $2}'`
L[4]=`grep "level5" ORBITAL_INPUT | awk -F "level5" '{print $2}'`
L[5]=`grep "level6" ORBITAL_INPUT | awk -F "level6" '{print $2}'`
L[6]=`grep "level7" ORBITAL_INPUT | awk -F "level7" '{print $2}'`
L[7]=`grep "level8" ORBITAL_INPUT | awk -F "level8" '{print $2}'`
L[8]=`grep "level9" ORBITAL_INPUT | awk -F "level9" '{print $2}'`

# (0.1.10)get some parameters for METROPOLIS
Start_tem_S_in=`grep "Start_tem_S" ORBITAL_INPUT | awk -F "Start_tem_S" '{print $0}' | awk '{print $2}'`
if ( test $Start_tem_S_in != " ") 
then
Start_tem_S=$Start_tem_S_in
else
Start_tem_S=1.0e-4                            #default
fi


Start_tem_K_in=`grep "Start_tem_K" ORBITAL_INPUT | awk -F "Start_tem_K" '{print $0}' | awk '{print $2}'`
if ( test $Start_tem_K_in != " " )
then
Start_tem_K=$Start_tem_K_in
else
Start_tem_K=1.0e-2                            #default
fi

Step_S_in=`grep "Step_S" ORBITAL_INPUT | awk -F "Step_S" '{print $0}' | awk '{print $2}'`
if ( test $Step_S_in != " " )
then
Step_S=$Step_S_in
else
Step_S=20                                     #default
fi

Step_K_in=`grep "Step_K" ORBITAL_INPUT | awk -F "Step_K" '{print $0}' | awk '{print $2}'`
if ( test $Step_K_in != " ")
then
Step_K=$Step_K_in
else
Step_K=15                                     #default
fi

Delta_kappa_in=`grep "Delta_kappa" ORBITAL_INPUT | awk -F "Delta_kappa" '{print $0}' | awk '{print $2}'`
Delta_kappa=$Delta_kappa_in
#echo "Delta_kappa=$Delta_kappa"

# (0.1.11) calculate the number of different dimers or trimers.
info=`grep "Dis" ORBITAL_INPUT | awk -F "Dis" '{print $2}'`
echo "info=$info"
dimers_number=`echo "$info" | awk '// {print NF}'`
echo "dimers_number=$dimers_number"

# (0.1.12) get the ecut
ecut=`grep "Ecut" ORBITAL_INPUT | awk -F "Ecut" '{print $0}' | awk '{print $2}'`

# (0.1.13) get the info about rcut,pseudo
info_r=`grep "Rcut" ORBITAL_INPUT | awk -F "Rcut" '{print $0}' | awk '{print $2}'`
rcut_number=`echo "$info_r" | awk '// {print NF}'`
echo "ecut=$ecut, info_r=$info_r, rcut_number=$rcut_number"

# (0.1.14) get the pseudopotential
pseudofile=`grep "Pseudo " ORBITAL_INPUT | awk -F "Pseudo " '{print $0}' | awk '{print $2}'`
echo "pseudo=$pseudofile"

# (0.1.15) get the smearing
smearing_sigma=`grep "smearing_sigma " ORBITAL_INPUT | awk -F "smearing_sigma " '{print $0}' | awk '{print $2}'`
echo "smearing_sigma=$smearing_sigma"

#-----------------------------------------------------------------
#
# (1) big cicle, cicle of targets
#
#-----------------------------------------------------------------

# (1.1) output which element you want to calculate

echo "--------------------------> $element"

# (1.2) make the dir, the name is 'name'
if ( test -d $name )
then 	 
	echo "The dir exist: $name"
else
	echo "make dir: $name"
	mkdir $name
fi

# (1.3) enter the name dir
cd $name

# (1.4) another big cicle come: the rcut cicle. 
count_r=1
while [ $count_r -le $rcut_number ]
do
	
	# (1.4.1)
	rcut=`echo "$info_r" | awk '{print $'$count_r'}' `
	echo "rcut=$rcut"


	# (1.4.2) enter the third big cicle : the dimer distance cicle.
	count=1
	while [ $count -le $dimers_number ]
	do

	# (1.4.2.0) calculate the distance of dimers
	dis=`echo "$info" | awk '{print $'$count'}' ` 
        dis1=$(echo "scale=5;$dis * 0.86603 "|bc)
        dis2=$(echo "scale=5;$dis * 0.5     "|bc)
	dis3=$(echo "scale=5;$dis * 0.81649 "|bc)
	dis4=$(echo "scale=5;$dis * 0.28867 "|bc)
        echo "dis=$dis"


# (1.4.2.1) get the" structures"
if ( test $element = Na -o $element = Li -o $element = K -o $element = Ca )
then
echo "use trimer"
na=3
cat > $name.stru << EOF
ATOMIC_SPECIES
$element $mass $pseudofile
LATTICE_CONSTANT
$lat0  // add lattice constant(a.u.)
LATTICE_VECTORS
1 0 0
0 1 0
0 0 1
ATOMIC_POSITIONS
Cartesian_angstrom  //Cartesian or Direct coordinate.
$element //Element Label
0.0     //starting magnetism
3       //number of atoms
0.0     0.0     0.0     0   0   0  // crystal coor.
0.0     0.0     $dis    0   0   0
0.0     $dis1   $dis2   0   0   0
EOF

#elif ( test $element = Mn -o $element = Fe ) 
#then
#echo "use tetramer"
#na=4
#cat > $name.stru << EOF
#ATOMIC_SPECIES
#$element $mass $pseudofile
#LATTICE_CONSTANT
#$lat0  // add lattice constant(a.u.)
#LATTICE_VECTORS
#1 0 0
#0 1 0
#0 0 1
#ATOMIC_POSITIONS
#Cartesian_angstrom  //Cartesian or Direct coordinate.
#$element //Element Label
#0.0     //starting magnetism
#4       //number of atoms
#0.0     0.0     0.0     0   0   0  // crystal coor.
#0.0     0.0     $dis    0   0   0 
#0.0     $dis1   $dis2   0   0   0 
#$dis3   $dis4   $dis2   0   0   0 
#EOF

else
echo "use dimer"
na=2
cat > $name.stru << EOF
ATOMIC_SPECIES
$element $mass $pseudofile
LATTICE_CONSTANT
$lat0  // add lattice constant(a.u.)
LATTICE_VECTORS
1 0 0
0 1 0
0 0 1
ATOMIC_POSITIONS
Cartesian_angstrom  //Cartesian or Direct coordinate.
$element //Element Label
0.0     //starting magnetism
2       //number of atoms
0.0     0.0     0.0     0   0   0  // crystal coor.
0.0     0.0     $dis    0   0   0
EOF
fi

# (1.4.2.2) get KPOINTS
cat > KPOINTS << EOF
K_POINTS
0
Gamma
1 1 1 0 0 0
EOF

# (1.4.2.3) get INPUTw
cat > INPUTw << EOF
WANNIER_PARAMETERS
rcut 10
out_spillage 2
EOF


# (1.4.2.4) get INPUTs
cat > INPUTs << EOF
INPUT_ORBITAL_INFORMATION
<SPHERICAL_BESSEL>
1           // smooth or not
0.1         // smearing_sigma
$ecut       // energy cutoff for spherical bessel functions(Ry)
$rcut       // cutoff of wavefunctions(a.u.)
1.0e-12     // tolerence
</SPHERICAL_BESSEL>
EOF


# (1.4.2.5) get INPUT
cat > INPUT << EOF
INPUT_PARAMETERS
suffix              $element-$rcut-$dis
latname             $element-$rcut-$dis
stru_file           $name.stru
pseudo_dir          $pseudo_dir
kpoint_file         KPOINTS
wannier_card        INPUTw
calculation         scf
ntype               1
nspin               1
lmaxmax             $maxL

symmetry            0
nbands             	$nbands 

ecutwfc             $ecut
scf_thr                 2.0e-8  // about iteration
scf_nmax               1000

smearing_method            gauss
smearing_sigma               $smearing_sigma

mixing_type         pulay       // about charge mixing
mixing_beta         0.4
mixing_ndim         8
printe				1
EOF

let count++

# (1.4.2.6)
#test -e ../node_openmpi && cp ../node_openmpi .

#-------------
#on Dirac
#-------------
#/opt/openmpi/bin/mpirun -np $cpu_number -machinefile node_openmpi $exe

#-------------
#on Einstein
#-------------
#mpiexec -np $cpu_num -machinefile node_openmpi $EXE_pw
#mpiexec -np $1 -machinefile $EXE_pw
mpirun -np $cpu_num $EXE_pw
#mpirun -np $cpu_num -hostfile $hostfpath  $EXE_pw

#$EXE_pw 

#mpiexec -n 12  -machinefile $PBS_NODEFILE $EXE_pw >> Log.txt
# end (1.4.2), dimer distace cicle
done

# (1.4.3) mkdir of rcut
test -d $rcut || mkdir $rcut
# (1.4.3.1)
cd $rcut

# (1.4.3.2) prepare for the INPUT file
cat > INPUT << EOF
<PW_QS>
1                       // if or not calculate the spillage. 1/0
0                       // restart or not. 1/0
1                       // if or not output the file. 1/0
$dimers_number          // number of structures.
EOF

# (1.4.3.3) input the file names
count_files=1
while [ $count_files -le $dimers_number ]
do
dis=`echo "$info" | awk '{print $'$count_files'}' ` 
cat >> INPUT << EOF
../$element-$rcut-$dis.$lat0.dat
EOF
let count_files++
done

cat >> INPUT << EOF
</PW_QS>

<PARALLEL>
1            //number of k points
1			// number of pools
</PARALLEL>

<BLOCK_NE>
100
</BLOCK_NE>

<SCHEME_VALUE>
2
<SCHEME_VALUE>

<METROPOLIS>
$Start_tem_S              // Start temparature for spillage
0.8                 // Cooling rate
$Step_S                   // Number of temperatures(spillage)
600                // Number of steps per temparature

$Start_tem_K              // start temperature for kinetic energy
0.8                 // Cooling rate
$Step_K                  // Number of temperatures(kinetical)
600                // Number of steps per temparature

$Delta_kappa        // Delta kappa
50                  // Selectly output information

100                 // Acceptance_steps
0.4                 // Acceptance_high
0.2                 // Acceptance_low

100                   // Max kinetic energy(Rydberg).
0.01                // 'dr' for kinetic minimized.
1                   // 1: Kin 2: Ecut
</METROPOLIS>


<BANDS>
1                   // to control the number of bands.(Yes1/No0)
1                   // int, the start band index(>0).
$ref_bands          // int, the ed band index(<provied bands).
</BANDS>
EOF

cat >> INPUT << EOF
<OPTIMIZE>
$Level                 // Number of levels.
label / na / skip / lmax / each L /
EOF

for((i=0;i<$Level;i++))
do
cat >> INPUT << EOF
$id $na new ${L[i]}
EOF
done

cat >> INPUT << EOF
</OPTIMIZE>
EOF

cat >> INPUT << EOF
<PLOT>
0.01    //dr(a.u.) of uniform mesh. Attention!!dr will affect kinetic energy minmized largely.
-6      //xmin
1       //zed, chosen as valence charge.
0.01    //dx
6.0     //xmax
</PLOT>

<CAL_C4>
0
2
./FILE/Si-S.ORBITAL
0
./FILE/Si-P.ORBITAL
1
</CAL_C4>

<TEST>
0                       // 'yes' to do this.
14.0                    // rcut, only useful for test program
0.01                    // dr, for simpson integral
2                       // test eigenvalue index
2                       // lmax
</TEST>

EOF

#mpiexec -n 1 -machinefile $PBS_NODEFILE $EXE_orbital >> Log.txt
$EXE_orbital 
#mpirun -np cpu_num $EXE_orbital

# (1.4.3.5)
cd ..

let count_r++
# end(1.4), rcut cicle
done

# (1.5) exit the name dir.
cd ..

done
# end targets cicle

unset OMP_NUM_THREADS
