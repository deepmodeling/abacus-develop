#!/bin/bash
# author mohan 
# edit by Pengfei Li 2013-6-4
# maintained by WenShuai Zhang (since 2016-11-30) with changelog:
# 1. fix bugs when reading "ORBITAL_INPUT" for lines begin with "#" or " ".
# 2. add support for generating orbital with multi different types of structures together, 
#    such as dimer STRU & tetramer STRU and so on.
# 3. add support for opt_orb_pytorch_dpsi generator

args=$@
if [ "$args" == "--help" -o "$args" == "-h" ] ; then 
    echo " Usage:    ./ThisScript.sh  <Input File Name> "
    echo " Example:  ./generate_orbital.sh  ORBITAL_INPUT_MixSTRU "
    echo " If no <Input File Name>, This script will use the default name for input file: ORBITAL_INPUT "
    echo " "
    exit
elif [ ! -z $1 ] && [ -f $1 ] ; then
    InputFile=$1
elif [ -f 'ORBITAL_INPUT' ] ; then 
    echo " Can't find specified input file, use default: ORBITAL_INPUT "
    InputFile="ORBITAL_INPUT"
else
    echo " Can't find specified input file or default file: ORBITAL_INPUT"
    echo " see help by '-h' "
    echo " "
    exit 
fi
if [ -z $InputFile ]; then 
    echo " Can't find input file"
    exit;
fi
#
echo " "
echo " ********************************************************* "
echo " *                                                       * "
echo " *          Start to Generate Orbital for LCAO           * "
echo " *                                                       * "
echo " ********************************************************* "

function completelyNewCalc(){
                echo " "
                echo " Completely New SIA Calculation ... "
                mkdir -p Old
                if ( test -f "ORBITAL_RESULTS.txt" ); then
                    echo " Move Old Orbital files: ORBITAL_RESULTS.txt to Old/"  
                    mv "ORBITAL_RESULTS.txt" Old/  
                    echo " Move Old INPUT   files: INPUT to Old/"
                    mv  "INPUT" Old/
                fi
}

time_start=`date +%s`
#-----------------------------------------------------------------
#
# (0.0) fixed  parameters 
#
#-----------------------------------------------------------------
mass=1                                         #not used yet.
lat0=20                                        #in a.u


echo "      InputFile = $InputFile"
#-----------------------------------------------------------------
# (0.0)   get EXE_mpi
#-----------------------------------------------------------------
EXE_mpi=`grep -E "^\s*EXE_mpi" $InputFile | awk -F "EXE_mpi" '{print $2}' `
EXE_mpi=`echo $EXE_mpi`
if [ -z "$EXE_mpi" ]; then
    cpu_num=8  
    hostfpath="myhosts"
    #hostfpath="/home/nic/wszhang/eclipse_project/mesia_dft/delta_test/delta_dft/cifs2deltaDFT/myhosts"
    grep slots $hostfpath > /dev/null
    if [ $? = 0 ] ; then        # has "slots" in file: $hostfpath
        cpu_num=`cat $hostfpath |cut -d "=" -f 2-3|awk '{sum += $1};END {print sum}'`
    else
        cpu_num=`cat $hostfpath |wc -l`
    fi
    echo "        cpu_num = $cpu_num"
    echo "      hostfpath = $hostfpath"
    EXE_mpi="mpirun -np $cpu_num -hostfile ../$hostfpath " 
fi
echo "        EXE_mpi = $EXE_mpi " 
Host_List=( `$EXE_mpi hostname 2>/dev/null` )
Host1_List=( `echo ${Host_List[@]} | grep -o "$Host_List"` )
Host1_NCore=${#Host1_List[@]} 
echo "    Host1_NCore = $Host1_NCore" 
if [ $Host1_NCore -lt 1 ]; then 
    echo " error at getting Host1_NCore, program exiting ... "
    exit 1
fi


#-----------------------------------------------------------------
# (0.1)   input parameters
#-----------------------------------------------------------------

# (0.1.1) get exe
EXE_pw=`grep -E "^\s*EXE_pw" $InputFile | awk -F "EXE_pw" '{print $0}' | awk '{print $2}'`
echo "         EXE_pw = $EXE_pw " 

# (0.1.2)get SIA
EXE_orbital=`grep -E "^\s*EXE_orbital" $InputFile | awk -F "EXE_orbital" '{print $0}' | awk '{print $2}'`
chmod +x "$EXE_orbital"


# (0.1.3)get the targets element and id
targets=`grep -E "^\s*targets" $InputFile | awk -F "targets" '{print $0}' | awk '{print $2}'`
echo "        targets = $targets" 



#-----------------------------------------------------------------
# (1) big cicle, cicle of targets
#-----------------------------------------------------------------
for name in $targets; 
do

#echo "name=$name" 
element=`echo "$name" | awk -F "_" '{print $2}'`
#echo "element=$element" 
id=`echo "$name" | awk -F "_" '{print $1}'`
#echo "id=$id" 

# (0.1.6)get the pseudo_dir
Pseudo_dir=`grep -E "^\s*Pseudo_dir" $InputFile | awk -F "Pseudo_dir" '{print $0}' | awk '{print $2}'`
echo "     Pseudo_dir = $Pseudo_dir " 

# (0.1.12) get the ecut
ecut=`grep -E "^\s*Ecut" $InputFile | awk -F "Ecut" '{print $0}' | awk '{print $2}'`
echo "           ecut = $ecut"

# (0.1.13) get the info about rcut,pseudo
#info_r=`grep -E "^\s*Rcut" $InputFile | awk -F "Rcut" '{print $0}' | awk '{print $2}'`
#rcut_number=`echo "$info_r" | awk '// {print NF}'`
info_r=( `grep -E "^\s*Rcut" $InputFile | awk -F "Rcut|#" '{print $2 }'` )
rcut_number=${#info_r[@]}
echo "    rcut_number = $rcut_number, info_r = ( ${info_r[@]} )"

# (0.1.14) get the pseudopotential
pseudofile=`grep -E "^\s*Pseudo " $InputFile | awk -F "Pseudo " '{print $0}' | awk '{print $2}'`
echo "         Pseudo = $pseudofile"

# (0.1.15) get the smearing
smearing_sigma=`grep -E "^\s*smearing_sigma " $InputFile | awk -F "smearing_sigma " '{print $0}' | awk '{print $2}'`
echo "          smearing_sigma = $smearing_sigma"


#
# (0.1.7)get maxL S P D
maxL=`grep -E "^\s*maxL" $InputFile | awk -F "maxL" '{print $0}' | awk '{print $2}'`

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


# (0.x.x) check info (include Level) for each STRU 
nSTRU=`grep -E "^\s*BLSTRU" $InputFile | wc -l`
#nSTRU=`grep -o "^\s*ListSTRU\s*[^#]*" W/ORBITAL_INPUT_DZP |wc -w |awk '{print $1-1}'`
if [[ "$nSTRU" == "0" ]] ; then 
    nSTRU=`grep -E "^\s*Dis" $InputFile | wc -l`
    if [ "$nSTRU" != "1" ]; then  
        echo "error: input file has more than one line for Dis argument, please use the BLSTRU* "
    fi
fi
echo "          nSTRU = $nSTRU"
#
LevelBegin[0]=0
EndLevel[0]=0
SkipSTRU[0]=0
RestartSTRU[0]=0
ListSTRU[0]=" "
for((iSTRU=1;iSTRU<=$nSTRU;iSTRU++))
do
    ListSTRU[iSTRU]=`grep -E "^\s*ListSTRU " $InputFile |awk -F "ListSTRU" '{print $0}' \\
                |awk -v iSTRU=$iSTRU '{print $(iSTRU+1) }'`
    ListSTRU[iSTRU]=${ListSTRU[iSTRU]:-dimer}
    echo "        STRU[$iSTRU] = ${ListSTRU[iSTRU]}"


    BLSTRU[iSTRU]=`grep -E "^\s*BLSTRU$iSTRU" $InputFile | awk -F "BLSTRU$iSTRU" '{print $2}'`
    if [[ -z "${BLSTRU[iSTRU]}" ]] && [[ "$iSTRU=1" ]] ; then 
        BLSTRU[iSTRU]=`grep -E "^\s*Dis" $InputFile | awk -F "Dis" '{print $2}'`
    fi
    BL_number[iSTRU]=`echo "${BLSTRU[iSTRU]}" | awk '// {print NF}'`
    # (0.1.11) calculate the number of different dimers or trimers.
    #BLSTRU=`grep "Dis1" $InputFile | awk -F "Dis1" '{print $2}'`
    #BL_number=`echo "$BLSTRU" | awk '// {print NF}'`
    echo "   BL_number[$iSTRU] = ${BL_number[iSTRU]}"
    echo "      BLSTRU[$iSTRU] ="  ${BLSTRU[iSTRU]} "  :Bond Length for each STRU"
    
    EndLevel[iSTRU]=`grep -E "^\s*Level" $InputFile |awk -F "Level" '{print $0}' \\
                |awk -v iSTRU=$iSTRU '{print $(iSTRU+1) }'`
    echo "    EndLevel[$iSTRU] = ${EndLevel[iSTRU]}"
    BeginLevel[iSTRU]=`grep -E "^\s*BeginLevel" $InputFile |awk -F "BeginLevel" '{print $0}' \\
                |awk -v iSTRU=$iSTRU '{print $(iSTRU+1) }'`
    echo "  BeginLevel[$iSTRU] = ${BeginLevel[iSTRU]:-auto/default}"
    
    # (0.1.4)get the nbands
    nbands[iSTRU]=`grep -E "^\s*nbands" $InputFile | awk -F "nbands" '{print $0}' \\
                    |awk -v iSTRU=$iSTRU '{print $(iSTRU+1) }'`
    echo "      nbands[$iSTRU] = ${nbands[iSTRU]}"
    
    # (0.1.5)get the ref_bands
    ref_bands[iSTRU]=`grep -E "^\s*ref_bands" $InputFile  | awk -F "$ref_bands" '{print $0}' \\
                    |awk -v iSTRU=$iSTRU '{print $(iSTRU+1) }'`
    echo "   ref_bands[$iSTRU] = ${ref_bands[iSTRU]}"

    RestartSTRU[iSTRU]=`grep -E "^\s*RestartSTRU" $InputFile \\
                    | awk -F "$RestartSTRU" '{print $0}' \\
                    | awk -v iSTRU=$iSTRU '{print $(iSTRU+1) }'`
    echo " RestartSTRU[$iSTRU] = ${RestartSTRU[iSTRU]:-auto/default}"
    #if [ ! -n "${RestartSTRU[iSTRU]}" ]; then
    #    RestartSTRU[iSTRU]=0
    #    echo " set RestartSTRU[$iSTRU]=0 "
    #fi


    SkipSTRU[iSTRU]=0
    if ( test -n "`grep -E "^\s*SkipSTRU" $InputFile`" ); then
        SkipSTRU[iSTRU]=`grep -E "^\s*SkipSTRU" $InputFile  \\
                    | awk -F "$SkipSTRU" '{print $0}' \\
                    | awk -v iSTRU=$iSTRU '{print $(iSTRU+1) }'`
    fi
    echo "    SkipSTRU[$iSTRU] = ${SkipSTRU[iSTRU]}"

done  # first cicle of iSTRU 
if [ "$nSTRU" == "1" ]; then 
    SkipSTRU[1]=0 
fi

# (0.1.8)get the level
#Level=`grep "Level" $InputFile | awk -F "level" '{print $0}' | awk '{print $2}'`
#echo "__Level=$Level"


# (0.1.9)get every level`s lmax s p d f g
Llevels[1]=`grep -E "^\s*level1" $InputFile | awk -F "level1" '{print $2}'`
Llevels[2]=`grep -E "^\s*level2" $InputFile | awk -F "level2" '{print $2}'`
Llevels[3]=`grep -E "^\s*level3" $InputFile | awk -F "level3" '{print $2}'`
Llevels[4]=`grep -E "^\s*level4" $InputFile | awk -F "level4" '{print $2}'`
Llevels[5]=`grep -E "^\s*level5" $InputFile | awk -F "level5" '{print $2}'`
Llevels[6]=`grep -E "^\s*level6" $InputFile | awk -F "level6" '{print $2}'`
Llevels[7]=`grep -E "^\s*level7" $InputFile | awk -F "level7" '{print $2}'`
Llevels[8]=`grep -E "^\s*level8" $InputFile | awk -F "level8" '{print $2}'`
Llevels[9]=`grep -E "^\s*level9" $InputFile | awk -F "level9" '{print $2}'`


# (0.1.10)get some parameters for METROPOLIS
Start_tem_S_in=`grep -E "^\s*Start_tem_S" $InputFile  \\
            | awk -F "Start_tem_S" '{print $0}' | awk '{print $2}'`
if ( test $Start_tem_S_in != " ") 
then
Start_tem_S=$Start_tem_S_in
else
Start_tem_S=1.0e-4                            #default
fi


Start_tem_K_in=`grep -E "^\s*Start_tem_K" $InputFile \\
            | awk -F "Start_tem_K" '{print $0}' | awk '{print $2}'`
if ( test $Start_tem_K_in != " " )
then
Start_tem_K=$Start_tem_K_in
else
Start_tem_K=1.0e-2                            #default
fi

Step_S_in=`grep -E "^\s*Step_S" $InputFile \\
        | awk -F "Step_S" '{print $0}' | awk '{print $2}'`
if ( test $Step_S_in != " " )
then
Step_S=$Step_S_in
else
Step_S=30                                     #default
fi

Step_K_in=`grep -E "^\s*Step_K" $InputFile \\
        | awk -F "Step_K" '{print $0}' | awk '{print $2}'`
if ( test $Step_K_in != " ")
then
Step_K=$Step_K_in
else
Step_K=20                                     #default
fi

Delta_kappa_in=`grep -E "^\s*Delta_kappa" $InputFile \\
            | awk -F "Delta_kappa" '{print $0}' | awk '{print $2}'`
Delta_kappa=$Delta_kappa_in
#echo "    Delta_kappa=$Delta_kappa"



# (1.1) output which element you want to calculate
echo " -------------------------------------------------------> $element"

# (1.2) make the dir, the name is 'name'
if ( test -d $name )
then 	 
	echo " The dir exist: $name"
else
	echo " Make dir: $name"
	mkdir $name
fi

# (1.3) enter the name dir
cd $name

## (1.4) another big cicle come: the rcut cicle. 
count_r=1
while [ $count_r -le $rcut_number ]
do
	### (1.4.1)
	#rcut=`echo "$info_r" | awk '{print $'$count_r'}' `
    rcut=${info_r[count_r-1]} 
	echo "    |run  cicle: rcut=$rcut"



    ### (1.4.1.5) enter the third big cicle: iSTRU  
    echo "    nSTRU = $nSTRU "
    for((iSTRU=1;iSTRU<=$nSTRU;iSTRU++))   
    do
        if ( test ${SkipSTRU[iSTRU]} -eq 1 ); then 
        echo "        |skip cicle: iSTRU=$iSTRU"
        continue;
        fi
	    echo "        |run  cicle: iSTRU=$iSTRU"



        ### (1.4.3) mkdir of rcut
        test -d $rcut || mkdir $rcut
        # (1.4.3.1)
        cd $rcut
        
        iSTRULeft=`expr $iSTRU \- 1`
        if [ "${BeginLevel[iSTRU]}" == "" ] ; then
            BeginLevel[$iSTRU]=$((${EndLevel[iSTRULeft]}+1))
            echo -e "\n not found BeginLevel[$iSTRU], use: EndLevel[iSTRULeft]+1"
        fi
        echo -e "\n BeginLevel[iSTRU]=${BeginLevel[iSTRU]},  EndLevel[iSTRU]=${EndLevel[iSTRU]} "
        


        ### set if restart from previous SIA runs 
        ### if ( test SkipSTRU[`expr $iSTRU - 1`] -eq 1 ) ; then
        #ifRestart=${RestartSTRU[$iSTRU]}
        if [ "${RestartSTRU[$iSTRU]}" == "" ] ; then
            if [ $iSTRU -gt 1 ] ; then
                RestartSTRU[$iSTRU]=1
            else
                RestartSTRU[$iSTRU]=0
            fi
        fi


        echo " RestartSTRU[$iSTRU] = ${RestartSTRU[$iSTRU]} "
        if [ ${RestartSTRU[$iSTRU]} -eq 0 ] ; then 
                completelyNewCalc
        else
            if [ $iSTRU -gt 1 ]; then

                echo " "
                echo " Current *.dat/*.txt ... will be considered previous calculation results of STRU${iSTRULeft} "
                #
                #if [ -f "ORBITAL_RESULTS.txt" ] ; then
                #    echo " Found file: ORBITAL_RESULTS.txt, continue ... "
                #else
                #    echo " Can't find: ORBITAL_RESULTS.txt, exiting ... "
                #    exit
                #fi
                # 
                # 
                if [ ! -f "STRU${iSTRULeft}.ORBITAL_RESULTS.txt" ]; then
                    echo " Move Previous Orbital files and Rename as STRU${iSTRULeft}.*"
                    #
                    if ( test -f "ORBITAL_RESULTS.txt" ); then
                        mv  "ORBITAL_RESULTS.txt"  "STRU${iSTRULeft}.ORBITAL_RESULTS.txt"
                    fi
                    #
                    if ( test -f "INPUT" ); then
                        mv  "INPUT" "STRU${iSTRULeft}.INPUT"
                    fi
                    #
                    if ( test -f "ORBITAL_${id}U.dat" ); then
                        mv  "ORBITAL_${id}U.dat" "STRU${iSTRULeft}.ORBITAL_${id}U.dat"
                    fi
                    #
                    if ( test -f "ORBITAL_${id}L.dat" ); then
                        mv  "ORBITAL_${id}L.dat" "STRU${iSTRULeft}.ORBITAL_${id}L.dat"
                    fi
                    #
                    if ( test -f "ORBITAL_ECUT.txt" ); then
                        mv  "ORBITAL_ECUT.txt"  "STRU${iSTRULeft}.ORBITAL_ECUT.txt"
                    fi
                    #
                    if ( test -f "ORBITAL_KINETIC.txt" ); then
                        mv  "ORBITAL_KINETIC.txt"  "STRU${iSTRULeft}.ORBITAL_KINETIC.txt"
                    fi
                    #
                    if ( test -f "ORBITAL_PLOTL.dat" ); then
                        mv  "ORBITAL_PLOTL.dat"  "STRU${iSTRULeft}.ORBITAL_PLOTL.dat"
                    fi
                    #
                    if ( test -f "ORBITAL_PLOTU.dat" ); then
                        mv  "ORBITAL_PLOTU.dat"  "STRU${iSTRULeft}.ORBITAL_PLOTU.dat"
                    fi
                    #
                    if ( test -f "ORBITAL_PLOTUK.dat" ); then
                        mv  "ORBITAL_PLOTUK.dat"  "STRU${iSTRULeft}.ORBITAL_PLOTUK.dat"
                    fi
                    #
                    if ( test -f "running_1.txt" ); then
                        mv  "running_1.txt" "STRU${iSTRULeft}.running_1.txt"
                    fi
                    #
                fi
                # 
                if [ -f "STRU${iSTRULeft}.ORBITAL_RESULTS.txt" ] ; then
                    echo " Found file: STRU${iSTRULeft}.ORBITAL_RESULTS.txt, copy as ORBITAL_RESULTS.txt ... " 
                    cp -ap "STRU${iSTRULeft}.ORBITAL_RESULTS.txt"  "ORBITAL_RESULTS.txt"  
                else
                    echo " Not found file: STRU${iSTRULeft}.ORBITAL_RESULTS.txt, exiting ... "
                    exit 
                fi
            else 
                echo " Current *.dat/*.txt ... will be considered previous calc. results of STRU${iSTRU} "
                echo " Before SIA Calculation: mv ... & cp ... "
                echo " Move: INPUT/*.dat/*.txt to Old/ "
                mkdir -p Old
                #
                if ( test -f "INPUT" ); then
                    mv  "INPUT"                    "Old/INPUT"
                fi
                if ( test -f "ORBITAL_${id}U.dat" ); then
                    mv  "ORBITAL_${id}U.dat"            "Old/ORBITAL_${id}U.dat"
                fi
                if ( test -f "ORBITAL_${id}L.dat" ); then
                    mv  "ORBITAL_${id}L.dat"            "Old/ORBITAL_${id}L.dat"
                fi
                if ( test -f "ORBITAL_ECUT.txt" ); then
                    mv  "ORBITAL_ECUT.txt"              "Old/ORBITAL_ECUT.txt"
                fi
                if ( test -f "ORBITAL_KINETIC.txt" ); then
                    mv  "ORBITAL_KINETIC.txt"           "Old/ORBITAL_KINETIC.txt"
                fi
                if ( test -f "ORBITAL_PLOTL.dat" ); then
                    mv  "ORBITAL_PLOTL.dat"             "Old/ORBITAL_PLOTL.dat"
                fi
                if ( test -f "ORBITAL_PLOTU.dat" ); then
                    mv  "ORBITAL_PLOTU.dat"             "Old/ORBITAL_PLOTU.dat"
                fi
                if ( test -f "ORBITAL_PLOTUK.dat" ); then
                    mv  "ORBITAL_PLOTUK.dat"            "Old/ORBITAL_PLOTUK.dat"
                fi
                if ( test -f "running_1.txt" ); then
                    mv  "running_1.txt"                 "Old/running_1.txt"
                fi
                #
                echo " Copy ORBITAL_RESULTS.txt to Old/ "
                if ( test -f "ORBITAL_RESULTS.txt" ); then
                    cp -avp  "ORBITAL_RESULTS.txt"      "Old/ORBITAL_RESULTS.txt"
                fi
                #
                echo " Start (SIA) Calculation:"
            fi
            echo " "
            echo " Restart from Previous Result: ORBITAL_RESULTS.txt "
            if [ -f "ORBITAL_RESULTS.txt" ] ; then
                echo " Found file: ORBITAL_RESULTS.txt " 
            else
                echo " Not found file: ORBITAL_RESULTS.txt, exiting ... "
                exit 
            fi
        fi  # which:: if [ ${RestartSTRU[$iSTRU]} -eq 0 ] ; then 
        echo " "

        

        ### (1.4.3.5) exit the rcut dir
        cd ..



	    #### (1.4.2) enter the forth big cicle : Bond Length.
	    count=1
	    while [ $count -le ${BL_number[iSTRU]} ]
	    do

	        # (1.4.2.0) calculate the Bond Length for iSTRU
	        BL=`echo "${BLSTRU[iSTRU]}" | awk '{print $'$count'}' ` 
            dis1=$(echo "scale=5;$BL * 0.86603 "|bc)
            dis2=$(echo "scale=5;$BL * 0.5     "|bc)
	        dis3=$(echo "scale=5;$BL * 0.81649 "|bc)
	        dis4=$(echo "scale=5;$BL * 0.28867 "|bc)
            echo "            |run  cicle: BL=$BL"

if [ "${ListSTRU[iSTRU]}" == "dimer" ]; then
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
0.0     0.0     $BL    0   0   0
EOF

elif [ "${ListSTRU[iSTRU]}" == "trimer" ]; then
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
0.0     0.0     $BL    0   0   0
0.0     $dis1   $dis2   0   0   0
EOF

elif [ "${ListSTRU[iSTRU]}" == "tetramer" ]; then
na=4
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
4       //number of atoms
0.0     0.0     0.0     0   0   0  // crystal coor.
0.0     0.0     $BL    0   0   0 
0.0     $dis1   $dis2   0   0   0 
$dis3   $dis4   $dis2   0   0   0 
EOF

fi 
echo " na=$na"


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
spillage_outdir OUT.$element-$rcut-$BL
EOF


# (1.4.2.4) get INPUTs
# cat > INPUTs << EOF
# INPUT_ORBITAL_INFORMATION
# <SPHERICAL_BESSEL>
# 1           // smooth or not
# 0.1         // smearing_sigma
# $ecut       // energy cutoff for spherical bessel functions(Ry)
# $rcut       // cutoff of wavefunctions(a.u.)
# 1.0e-12     // tolerence
# </SPHERICAL_BESSEL>
# EOF


# (1.4.2.5) get INPUT
cat > INPUT << EOF
INPUT_PARAMETERS
suffix              $element-$rcut-$BL
stru_file           $name.stru
pseudo_dir          $Pseudo_dir
kpoint_file         KPOINTS
wannier_card        INPUTw
calculation         scf
ntype               1
nspin               1
lmaxmax             $maxL
bessel_nao_rcut     $rcut

symmetry            0
nbands              ${nbands[iSTRU]} 

ecutwfc             $ecut
scf_thr             1e-7  // about iteration
scf_nmax            1500

smearing_method     gauss
smearing_sigma      $smearing_sigma

mixing_type         pulay       // about charge mixing
mixing_beta         0.4
mixing_ndim         8
printe              1
EOF

let count++

        echo " pwd:"
        pwd 
        if [ ${RestartSTRU[$iSTRU]} -eq 2 ]; then #  grep -E "^\s*RestartSTRU" ../$InputFile  > /dev/null 2>&1 ; then 
            echo " Skip_Calculation: $EXE_mpi $EXE_pw"
        else 
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
            
            #exit
            #mpiexec -n 12  -machinefile $PBS_NODEFILE $EXE_pw >> Log.txt
            #mpirun -np $cpu_num $EXE_pw
            #mpirun -hostfile "../$hostfpath"  $EXE_pw
            #mpirun -np $cpu_num -hostfile "../$hostfpath"  $EXE_pw

            echo " $EXE_mpi $EXE_pw "
            $EXE_mpi $EXE_pw 
            
            #mv test.0.dat  OUT.$element-$rcut-$BL/orb_matrix.0.dat
            #mv test.1.dat  OUT.$element-$rcut-$BL/orb_matrix.1.dat
            #echo " Skip_Calculation: $EXE_mpi $EXE_pw"
        fi
        echo ""



        #### end cicle (1.4.2): Bond Length 
        done



        # (1.4.3.1)
        cd $rcut



        ###if [ "${EXE_orbital##*.}" != "py" ]; 
        if [ "${EXE_orbital:0-3:3}" != ".py" ]; 
        then 
            echo -e " Using Simulated Annealing Method \n" 

    ### (1.4.3.2) prepare for the INPUT file
cat > INPUT << EOF
<PW_QS>
1                       // if or not calculate the spillage. 1/0
$(( ${RestartSTRU[$iSTRU]} > 0 ))                       // restart or not. 1/0
1                       // if or not output the file. 1/0
${BL_number[iSTRU]}          // number of structures.
EOF

# (1.4.3.3) input the file names
count_files=1
while [ $count_files -le ${BL_number[iSTRU]} ]
do
    BL=`echo "${BLSTRU[iSTRU]}" | awk '{print $'$count_files'}' ` 
    cat >> INPUT << EOF
../OUT.$element-$rcut-$BL/orb_matrix.0.dat
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
${ref_bands[$iSTRU]}          // int, the ed band index(<provied bands).
</BANDS>

EOF

cat >> INPUT << EOF
<OPTIMIZE>
${EndLevel[$iSTRU]}             // Number of levels.
label / na / skip / lmax / each L /
EOF


for((i=1;i<=${EndLevel[iSTRU]};i++))
do
    if [ $i -ge ${BeginLevel[iSTRU]} -a $i -le ${EndLevel[iSTRU]} ] 
    then
        leveltype="new "
    else
        leveltype="skip"
    fi
#echo "leveltype=$leveltype"
cat >> INPUT << EOF
$id $na $leveltype ${Llevels[i]}
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














        elif [ "${EXE_orbital:0-12:12}" == "dpsi/main.py"  ]; then  # begin from if [ "${EXE_orbital:0-3:3}" != ".py" ]; 
            echo -e " Using PyTorch Gradient Method V2 (include dpsi) \n" 

            ### len(dis[info["input"]["element"]]),
            ### (1.4.3.2) prepare INPUT file in json for PyTorch program

cat > INPUT << EOF
{
    "file_list":
    {
        "origin":
        [  
EOF

count_files=1
while [ $count_files -lt ${BL_number[iSTRU]} ]
do
    BL=`echo "${BLSTRU[iSTRU]}" | awk '{print $'$count_files'}' ` 
    cat >> INPUT << EOF
            "../OUT.$element-$rcut-$BL/orb_matrix.0.dat",
EOF
    let count_files++
done
BL=`echo "${BLSTRU[iSTRU]}" | awk '{print $'$count_files'}' ` 
cat >> INPUT << EOF
            "../OUT.$element-$rcut-$BL/orb_matrix.0.dat"
        ],
        "linear":
        [
            [ 
EOF

count_files=1
while [ $count_files -lt ${BL_number[iSTRU]} ]
do
    BL=`echo "${BLSTRU[iSTRU]}" | awk '{print $'$count_files'}' ` 
    cat >> INPUT << EOF
            "../OUT.$element-$rcut-$BL/orb_matrix.1.dat",
EOF
    let count_files++
done
BL=`echo "${BLSTRU[iSTRU]}" | awk '{print $'$count_files'}' ` 
cat >> INPUT << EOF
            "../OUT.$element-$rcut-$BL/orb_matrix.1.dat"
            ]
        ]
EOF

cat >> INPUT << EOF

    },
    "info":
    {
        "Nt_all": [
            "$element"
            ],
        "Nu": {
            "$element": [
EOF


if [ ${RestartSTRU[$iSTRU]} -ge 1 ] ; then 
    lr_value=0.0001
else
    lr_value=0.01
fi
echo  " iSTRULeft: $iSTRULeft, iSTRU: $iSTRU "
echo  " BeginLevel[STRUs]: (${BeginLevel[@]}) "
echo  " EndLevel[STRUs]: (${EndLevel[@]}) "
#if   [ $ifRestart -eq 0 ] ; then
#    C_init_from_file="false"
#elif [ $ifRestart -eq 1 ] ; then
#    C_init_from_file="true"
#fi


LValueMax=0
for LValue in {0..4} ;
do 
    numL[$LValue]=0
done

for((i=1;i<=${EndLevel[iSTRU]};i++))
do
    #if [ $i -le ${EndLevel[iSTRULeft]} ]; 
    #then
    #    C_init_from_file="true"
    #    #echo " Level:$i, C_init_from_file = " $C_init_from_file
    #fi

    Llevels_i=( ${Llevels[i]} )
    for LValue in {0..4} ;
    do 
        nAdd=${Llevels_i[$LValue+1]}
        #echo " nAdd = $nAdd"
        if [ "$nAdd" != "" ]; then 
            numL[$LValue]=$(( ${numL[$LValue]} + $nAdd ))

            if [ $LValueMax -lt $LValue ]; then 
                LValueMax=$LValue
            fi
        fi
    done 
done
echo    " numL = ${numL[@]}, LValueMax = $LValueMax "

for((LValue=0; LValue<${LValueMax}; LValue++))
do
    cat >> INPUT << EOF
                ${numL[$LValue]}, 
EOF
done

    cat >> INPUT << EOF
                ${numL[$LValueMax]}
                ]
        },
        "Rcut": {
            "$element": $rcut 
        },
        "dr": {
            "$element": 0.01
        },
        "Ecut": {
            "$element": $ecut
        },
    	"lr":  $lr_value,
        "cal_T":    true,
        "cal_smooth": true
    },
	"weight":
	{
        "stru": [
EOF

for((i_BL=1; i_BL<${BL_number[iSTRU]}; i_BL++))
do
    cat >> INPUT << EOF
            1, 
EOF
done
    cat >> INPUT << EOF
            1
         ],
EOF


if [ 2 -ge 1 ] ; then     
    cat >> INPUT << EOF
        "bands_range": [
EOF
for((i_BL=1; i_BL<${BL_number[iSTRU]}; i_BL++))
do
    cat >> INPUT << EOF
            `printf %.0f ${ref_bands[iSTRU]}`, 
EOF
done
    cat >> INPUT << EOF
            `printf %.0f ${ref_bands[iSTRU]}`
         ]
EOF
else
    cat >> INPUT << EOF
        "bands_file": [
         ]
EOF
fi    


    cat >> INPUT << EOF
    },
    "C_init_info":
    {
EOF
C_init_file="ORBITAL_RESULTS.txt"
if [ ${RestartSTRU[$iSTRU]} -ge 1 ] ; then  ## [ "$C_init_from_file" == "true" ]; then
    cat >> INPUT << EOF
        "init_from_file": true,
        "C_init_file"   : "$C_init_file",
        "opt_C_read"    : false
EOF
    echo -e " init_from_file : true, \n C_init_file : $C_init_file, \n opt_C_read: false\n"
elif [ ${RestartSTRU[$iSTRU]} -eq 0 ] ; then
    cat >> INPUT << EOF
        "init_from_file": false
EOF
    echo -e " init_from_file : false \n "
else 
    exit
fi

    cat >> INPUT << EOF
    },
    "V_info": 
    {
        "same_band":		true,
        "init_from_file":	true
    }
}
EOF


export OMP_NUM_THREADS=$Host1_NCore
echo " Set OMP_NUM_THREADS = $OMP_NUM_THREADS "
echo " Python2 Version: " `which python2`
echo " Python3 Version: " `which python3`
echo "" 




        
        
        







        else # begin from: if [ "${EXE_orbital:0-3:3}" != ".py" ];
            echo -e " Using Previous PyTorch Gradient Method \n" 

            ### len(dis[info["input"]["element"]]),
            ### (1.4.3.2) prepare INPUT file in json for PyTorch program

cat > INPUT << EOF
{
    "file_list": [
EOF

count_files=1
while [ $count_files -lt ${BL_number[iSTRU]} ]
do
    BL=`echo "${BLSTRU[iSTRU]}" | awk '{print $'$count_files'}' ` 
    cat >> INPUT << EOF
        "../OUT.$element-$rcut-$BL/orb_matrix.0.dat",
EOF
    let count_files++
done
    BL=`echo "${BLSTRU[iSTRU]}" | awk '{print $'$count_files'}' ` 
cat >> INPUT << EOF
        "../OUT.$element-$rcut-$BL/orb_matrix.0.dat"
    ],
    "info": {
        "Nt_all": [
            "$element"
            ],
        "Nu": {
            "$element": [
EOF


echo  " iSTRULeft: $iSTRULeft, iSTRU: $iSTRU "
echo  " BeginLevel[STRUs]: (${BeginLevel[@]}) "
echo  " EndLevel[STRUs]: (${EndLevel[@]}) "
#if   [ $ifRestart -eq 0 ] ; then
#    C_init_from_file="false"
#elif [ $ifRestart -eq 1 ] ; then
#    C_init_from_file="true"
#fi

LValueMax=0
for LValue in {0..4} ;
do 
    numL[$LValue]=0
done

for((i=1;i<=${EndLevel[iSTRU]};i++))
do
    #if [ $i -le ${EndLevel[iSTRULeft]} ]; 
    #then
    #    C_init_from_file="true"
    #    #echo " Level:$i, C_init_from_file = " $C_init_from_file
    #fi

    Llevels_i=( ${Llevels[i]} )
    for LValue in {0..4} ;
    do 
        nAdd=${Llevels_i[$LValue+1]}
        #echo " nAdd = $nAdd"
        if [ "$nAdd" != "" ]; then 
            numL[$LValue]=$(( ${numL[$LValue]} + $nAdd ))

            if [ $LValueMax -lt $LValue ]; then 
                LValueMax=$LValue
            fi
        fi
    done 
done
echo    " numL = ${numL[@]}, LValueMax = $LValueMax "


for((LValue=0; LValue<${LValueMax}; LValue++))
do
    cat >> INPUT << EOF
                ${numL[$LValue]}, 
EOF
done
    cat >> INPUT << EOF
                ${numL[$LValueMax]}
                ]
        },
        "Nb_true": [
EOF

for((i_BL=1; i_BL<${BL_number[iSTRU]}; i_BL++))
do
    cat >> INPUT << EOF
            ${ref_bands[iSTRU]}, 
EOF
done
    cat >> INPUT << EOF
            ${ref_bands[iSTRU]}
            ],
        "weight": [
EOF

for((i_BL=1; i_BL<${BL_number[iSTRU]}; i_BL++))
do
    cat >> INPUT << EOF
            1, 
EOF
done
if [ ${RestartSTRU[$iSTRU]} -ge 1 ] ; then 
    lr_value=0.0001
else
    lr_value=0.01
fi
    cat >> INPUT << EOF
            1
            ],
        "Rcut": {
            "$element": $rcut 
        },
        "dr": {
            "$element": 0.01
        },
        "Ecut": {
            "$element": $ecut
        },
    	"lr":  $lr_value,
        "cal_T":    false,
        "cal_smooth": true
    },
    "C_init_info": {
EOF

C_init_file="ORBITAL_RESULTS.txt"
if [ ${RestartSTRU[$iSTRU]} -ge 1 ] ; then  ## [ "$C_init_from_file" == "true" ]; then
    cat >> INPUT << EOF
        "init_from_file": true,
        "C_init_file"   : "$C_init_file"
EOF
    echo -e " init_from_file : true \n C_init_file : $C_init_file \n "
elif [ ${RestartSTRU[$iSTRU]} -eq 0 ] ; then
    cat >> INPUT << EOF
        "init_from_file": false
EOF
    echo -e " init_from_file : false \n "
else 
    exit
fi

    cat >> INPUT << EOF
    },
    "V_info": {
        "same_band": true,
        "init_from_file": false
    }
}
EOF



export OMP_NUM_THREADS=$Host1_NCore
echo " Set OMP_NUM_THREADS = $OMP_NUM_THREADS "
echo " Python2 Version: " `which python2`
echo " Python3 Version: " `which python3`
echo "" 

        fi  # begin from: if [ "${EXE_orbital:0-3:3}" != ".py" ];
        #cat INPUT



        #mpiexec -n 1 -machinefile $PBS_NODEFILE $EXE_orbital >> Log.txt
        echo " Run $EXE_orbital"
        echo ""
        $EXE_orbital 
        #mpirun -np cpu_num $EXE_orbital



        if [ "${EXE_orbital:0-3:3}" == ".py" ]; then
            echo "" 
            unset  OMP_NUM_THREADS
            echo " unset OMP_NUM_THREADS "
            #... 
            #echo " Back to Python2: " `which python2`
            #echo " Back to Python3: " `which python3`
            sleep 1
        fi

        #exit

    ### (1.4.3.5) exit the rcut dir
    cd ..
    ### end cicle (1.4.1.5): iSTRU 
    done  



let count_r++
## end cicle (1.4): rcut
done



# (1.5) exit the name dir.
cd ..
done
# end cicle (1): targets

time_end=`date +%s`
time_passed=$(($time_end - $time_start)) 
echo -e " Total Time:   $time_passed  \n"

