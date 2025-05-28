#!/bin/bash

# mohan add 2025-05-03
# this compare script is used in different integrate tests
COMPARE_SCRIPT="../../integrate/tools/CompareFile.py"
SUM_CUBE_EXE="../../integrate/tools/sum_cube.exe"


sum_file(){
	line=`grep -vc '^$' $1`
	inc=1
	if ! test -z $2; then
		inc=$2
	fi	
	sum=0.0
	for (( num=1 ; num<=$line ; num+=$inc ));do
		value_line=(` sed -n "$num p" $1 | head -n1 `)
		colume=`echo ${#value_line[@]}`
		for (( col=0 ; col<$colume ; col++ ));do
			value=`echo ${value_line[$col]}`
			sum=`awk 'BEGIN {x='$sum';y='$value';printf "%.6f\n",x+sqrt(y*y)}'` 
		done
	done
	echo $sum
}


get_input_key_value(){
	key=$1
	inputf=$2
	value=$(awk -v key=$key '{if($1==key) a=$2} END {print a}' $inputf)
	echo $value
}


file=$1
#echo $1

# the command will ignore lines starting with #
calculation=`grep calculation INPUT | grep -v '^#' | awk '{print $2}' | sed s/[[:space:]]//g`

running_path=`echo "OUT.autotest/running_$calculation"".log"`
#echo $running_path

natom=`grep -En '(^|[[:space:]])TOTAL ATOM NUMBER($|[[:space:]])' $running_path | tail -1 | awk '{print $6}'`
has_force=$(get_input_key_value "cal_force" "INPUT")
has_stress=$(get_input_key_value "cal_stress" "INPUT")
has_dftu=$(get_input_key_value "dft_plus_u" "INPUT")
has_band=$(get_input_key_value "out_band" "INPUT")
has_dos=$(get_input_key_value "out_dos" "INPUT")
has_cond=$(get_input_key_value "cal_cond" "INPUT")
has_hs=$(get_input_key_value "out_mat_hs" "INPUT")
has_hs2=$(get_input_key_value "out_mat_hs2" "INPUT")
has_xc=$(get_input_key_value "out_mat_xc" "INPUT")
has_xc2=$(get_input_key_value "out_mat_xc2" "INPUT")
has_eband_separate=$(get_input_key_value "out_eband_terms" "INPUT")
has_r=$(get_input_key_value "out_mat_r" "INPUT")
has_lowf=$(get_input_key_value "out_wfc_lcao" "INPUT")
out_app_flag=$(get_input_key_value "out_app_flag" "INPUT")
has_wfc_r=$(get_input_key_value "out_wfc_r" "INPUT")
has_wfc_pw=$(get_input_key_value "out_wfc_pw" "INPUT")
out_dm=$(get_input_key_value "out_dm" "INPUT")
out_mul=$(get_input_key_value "out_mul" "INPUT")
gamma_only=$(get_input_key_value "gamma_only" "INPUT")
imp_sol=$(get_input_key_value "imp_sol" "INPUT")
run_rpa=$(get_input_key_value "rpa" "INPUT")
out_pot=$(get_input_key_value "out_pot" "INPUT")
out_elf=$(get_input_key_value "out_elf" "INPUT")
out_dm1=$(get_input_key_value "out_dm1" "INPUT")
get_s=$(get_input_key_value "calculation" "INPUT")
out_pband=$(get_input_key_value "out_proj_band" "INPUT")
toW90=$(get_input_key_value "towannier90" "INPUT")
has_mat_r=$(get_input_key_value "out_mat_r" "INPUT")
has_mat_t=$(get_input_key_value "out_mat_t" "INPUT") 
has_mat_dh=$(get_input_key_value "out_mat_dh" "INPUT")
has_scan=$(get_input_key_value "dft_functional" "INPUT")
out_chg=$(get_input_key_value "out_chg" "INPUT") 
has_ldos=$(get_input_key_value "out_ldos" "INPUT")
esolver_type=$(get_input_key_value "esolver_type" "INPUT")
rdmft=$(get_input_key_value "rdmft" "INPUT")
base=$(get_input_key_value "basis_type" "INPUT")
word_total_time="atomic_world"
symmetry=$(get_input_key_value "symmetry" "INPUT")
out_current=$(get_input_key_value "out_current" "INPUT")
test -e $1 && rm $1

#------------------------------------------------------------
# if NOT non-self-consistent calculations or linear response
#------------------------------------------------------------
is_lr=0
if [ ! -z $esolver_type ] && ([ $esolver_type == "lr" ] || [ $esolver_type == "ks-lr" ]); then
	is_lr=1
fi

#----------------------------
# total energy information
#----------------------------
if [ $calculation != "get_wf" ]\
&& [ $calculation != "get_pchg" ] && [ $calculation != "get_S" ]\
&& [ $is_lr == 0 ]; then
	etot=$(grep "ETOT_" "$running_path" | tail -1 | awk '{print $2}')
    #echo "etot = $etot"
	etotperatom=`awk 'BEGIN {x='$etot';y='$natom';printf "%.10f\n",x/y}'`
    #echo "etotperatom = $etotperatom"
    # put the results in file
	echo "etotref $etot" >>$1
	echo "etotperatomref $etotperatom" >>$1
fi

#----------------------------
# force information
# echo "hasforce:"$has_force
#----------------------------
if ! test -z "$has_force" && [ $has_force == 1 ]; then
	nn3=`echo "$natom + 3" |bc`
    # echo "nn3=$nn3"
    # check the last step result
    grep -A$nn3 "TOTAL-FORCE" $running_path |awk 'NF==4{print $2,$3,$4}' | tail -$natom > force.txt
	total_force=`sum_file force.txt`
    rm force.txt
	echo "totalforceref $total_force" >>$1
fi

#-------------------------------
# stress information
# echo "has_stress:"$has_stress
#-------------------------------
if ! test -z "$has_stress" && [  $has_stress == 1 ]; then
    grep -A6 "TOTAL-STRESS" $running_path| awk 'NF==3' | tail -3> stress.txt
	total_stress=`sum_file stress.txt`
	rm stress.txt
	echo "totalstressref $total_stress" >>$1
fi


#-------------------------------
# DOS information
# echo $total_charge
#-------------------------------
if ! test -z "$has_dos"  && [  $has_dos == 1 ]; then
	total_dos=`cat OUT.autotest/DOS1_smear.dat | awk 'END {print}' | awk '{print $3}'`
	echo "totaldosref $total_dos" >> $1
fi

#-------------------------------
# Onsager coefficiency
#-------------------------------
if ! test -z "$has_cond"  && [  $has_cond == 1 ]; then
	onref=refOnsager.txt
	oncal=OUT.autotest/Onsager.txt
	python3 $COMPARE_SCRIPT $onref $oncal 3 -com_type 0
    echo "CompareH_Failed $?" >>$1
	rm -f je-je.txt Chebycoef
fi

#-------------------------------
# echo $out_dm1
#-------------------------------
if ! test -z "$out_dm1"  && [  $out_dm1 == 1 ]; then
	dm1ref=dmrs1_nao.csr.ref
	dm1cal=OUT.autotest/dmrs1_nao.csr
	python3 $COMPARE_SCRIPT $dm1ref $dm1cal 8
	echo "CompareDM1_pass $?" >>$1
fi

#-------------------------------
# echo $out_pot1
#-------------------------------
if ! test -z "$out_pot"  && [  $out_pot == 1 ]; then
	pot1ref=pots1.cube.ref
	pot1cal=OUT.autotest/pots1.cube
	python3 $COMPARE_SCRIPT $pot1ref $pot1cal 3
	echo "ComparePot1_pass $?" >>$1
fi

#-------------------------------
#echo $out_pot2
#-------------------------------
if ! test -z "$out_pot"  && [  $out_pot == 2 ]; then
	pot1ref=pot_es.cube.ref
	pot1cal=OUT.autotest/pot_es.cube
	python3 $COMPARE_SCRIPT $pot1ref $pot1cal 8
	echo "ComparePot_pass $?" >>$1
fi

#-------------------------------
# Electron localized function
# echo $out_elf
#-------------------------------
if ! test -z "$out_elf"  && [  $out_elf == 1 ]; then
	elf1ref=refELF.cube
	elf1cal=OUT.autotest/ELF.cube
	python3 $COMPARE_SCRIPT $elf1ref $elf1cal 3
	echo "ComparePot1_pass $?" >>$1
fi

#-------------------------------
# Overlap matrix
# echo $get_s
#-------------------------------
if ! test -z "$get_s"  && [  $get_s == "get_S" ]; then
	sref=refSR.csr
	scal=OUT.autotest/SR.csr
	python3 $COMPARE_SCRIPT $sref $scal 8
	echo "CompareS_pass $?" >>$1
fi

#-------------------------------
# Partial band structure
# echo $out_pband
#-------------------------------
if ! test -z "$out_pband"  && [  $out_pband == 1 ]; then
	#pbandref=refPBANDS_1
	#pbandcal=OUT.autotest/PBANDS_1
	#python3 $COMPARE_SCRIPT $pbandref $pbandcal 8
	#echo "CompareProjBand_pass $?" >>$1
	orbref=refOrbital
	orbcal=OUT.autotest/Orbital
	python3 $COMPARE_SCRIPT $orbref $orbcal 8
	echo "CompareOrb_pass $?" >>$1
fi

#-------------------------------
# Wannier90 information
# echo $toW90
#-------------------------------
if ! test -z "$toW90"  && [  $toW90 == 1 ]; then
	amnref=diamond.amn
	amncal=OUT.autotest/diamond.amn
	mmnref=diamond.mmn
	mmncal=OUT.autotest/diamond.mmn
	eigref=diamond.eig
	eigcal=OUT.autotest/diamond.eig
	sed -i '1d' $amncal
	sed -i '1d' $mmncal
	python3 $COMPARE_SCRIPT $amnref $amncal 1 -abs 8
	echo "CompareAMN_pass $?" >>$1
	python3 $COMPARE_SCRIPT $mmnref $mmncal 1 -abs 8
	echo "CompareMMN_pass $?" >>$1
	python3 $COMPARE_SCRIPT $eigref $eigcal 8
	echo "CompareEIG_pass $?" >>$1
fi

#-------------------------------
# Total DOS
# echo total_dos
# echo $has_band
#-------------------------------
if ! test -z "$has_band"  && [  $has_band == 1 ]; then
	bandref=refBANDS_1.dat
	bandcal=OUT.autotest/BANDS_1.dat
	python3 $COMPARE_SCRIPT $bandref $bandcal 8
	echo "CompareBand_pass $?" >>$1
fi


#--------------------------------
# Hamiltonian and overlap matrix
# echo $has_hs
#--------------------------------
if ! test -z "$has_hs"  && [  $has_hs == 1 ]; then
	if ! test -z "$gamma_only"  && [ $gamma_only == 1 ]; then
                href=data-0-H.ref
                hcal=OUT.autotest/data-0-H
                sref=data-0-S.ref
                scal=OUT.autotest/data-0-S
        else
                href=data-1-H.ref
                hcal=OUT.autotest/data-1-H
                sref=data-1-S.ref
                scal=OUT.autotest/data-1-S
        fi

        python3 $COMPARE_SCRIPT $href $hcal 6
    echo "CompareH_pass $?" >>$1
    python3 $COMPARE_SCRIPT $sref $scal 8
    echo "CompareS_pass $?" >>$1
fi

#--------------------------------
# exchange-correlation potential 
#--------------------------------
if ! test -z "$has_xc"  && [  $has_xc == 1 ]; then
	if ! test -z "$gamma_only"  && [ $gamma_only == 1 ]; then
			xcref=k-0-Vxc.ref
			xccal=OUT.autotest/k-0-Vxc
	else
			xcref=k-1-Vxc.ref
			xccal=OUT.autotest/k-1-Vxc
	fi
	oeref=vxc_out.ref
	oecal=OUT.autotest/vxc_out.dat
	python3 $COMPARE_SCRIPT $xcref $xccal 4
	echo "CompareVXC_pass $?" >>$1
	python3 $COMPARE_SCRIPT $oeref $oecal 5
    echo "CompareOrbXC_pass $?" >>$1
fi

#--------------------------------
# exchange-correlation potential 
#--------------------------------
if ! test -z "$has_xc2"  && [  $has_xc2 == 1 ]; then
	xc2ref=Vxc_R_spin0.ref
	xc2cal=OUT.autotest/Vxc_R_spin0.csr
	python3 $COMPARE_SCRIPT $xc2ref $xc2cal 8
	echo "CompareVXC_R_pass $?" >>$1
fi

#--------------------------------
# separate terms in band enegy 
#--------------------------------
if ! test -z "$has_eband_separate"  && [  $has_eband_separate == 1 ]; then
	ekref=kinetic_out.ref
	ekcal=OUT.autotest/kinetic_out.dat
	python3 $COMPARE_SCRIPT $ekref $ekcal 4
	echo "CompareOrbKinetic_pass $?" >>$1
	vlref=vpp_local_out.ref
	vlcal=OUT.autotest/vpp_local_out.dat
	python3 $COMPARE_SCRIPT $vlref $vlcal 4
	echo "CompareOrbVL_pass $?" >>$1
	vnlref=vpp_nonlocal_out.ref
	vnlcal=OUT.autotest/vpp_nonlocal_out.dat
	python3 $COMPARE_SCRIPT $vnlref $vnlcal 4
	echo "CompareOrbVNL_pass $?" >>$1
	vhref=vhartree_out.ref
	vhcal=OUT.autotest/vhartree_out.dat
	python3 $COMPARE_SCRIPT $vhref $vhcal 4
	echo "CompareOrbVHartree_pass $?" >>$1
fi

#-----------------------------------
# Hamiltonian and overlap matrices
#-----------------------------------
#echo $has_hs2
if ! test -z "$has_hs2"  && [  $has_hs2 == 1 ]; then
    #python3 $COMPARE_SCRIPT data-HR-sparse_SPIN0.csr.ref OUT.autotest/data-HR-sparse_SPIN0.csr 8
    #echo "CompareHR_pass $?" >>$1
    python3 $COMPARE_SCRIPT data-SR-sparse_SPIN0.csr.ref OUT.autotest/data-SR-sparse_SPIN0.csr 8
    echo "CompareSR_pass $?" >>$1
fi

#-----------------------------------
#  <psi_i0 | r | psi_jR> matrix
#-----------------------------------
#echo $has_mat_r
if ! test -z "$has_mat_r"  && [  $has_mat_r == 1 ]; then
    python3 $COMPARE_SCRIPT data-rR-sparse.csr.ref OUT.autotest/data-rR-sparse.csr 8
    echo "ComparerR_pass $?" >>$1
fi

#-----------------------------------
#  <psi_i0 | T | psi_jR> matrix
#-----------------------------------
#echo $has_mat_t
if ! test -z "$has_mat_t"  && [  $has_mat_t == 1 ]; then
    python3 $COMPARE_SCRIPT data-TR-sparse_SPIN0.csr.ref OUT.autotest/data-TR-sparse_SPIN0.csr 8
    echo "ComparerTR_pass $?" >>$1
fi

#-----------------------------------
#  <psi_i0 | dH | psi_jR> matrix
#-----------------------------------
#echo $has_mat_dh
if ! test -z "$has_mat_dh"  && [  $has_mat_dh == 1 ]; then
    python3 $COMPARE_SCRIPT data-dHRx-sparse_SPIN0.csr.ref OUT.autotest/data-dHRx-sparse_SPIN0.csr 8
    echo "ComparerdHRx_pass $?" >>$1
    python3 $COMPARE_SCRIPT data-dHRy-sparse_SPIN0.csr.ref OUT.autotest/data-dHRy-sparse_SPIN0.csr 8
    echo "ComparerdHRy_pass $?" >>$1
    python3 $COMPARE_SCRIPT data-dHRz-sparse_SPIN0.csr.ref OUT.autotest/data-dHRz-sparse_SPIN0.csr 8
    echo "ComparerdHRz_pass $?" >>$1
fi

#---------------------------------------
# SCAN exchange-correlation information
#echo $has_scan
#---------------------------------------
if ! test -z "$has_scan"  && [  $has_scan == "scan" ] && \
       ! test -z "$out_chg" && [ $out_chg == 1 ]; then
    python3 $COMPARE_SCRIPT SPIN1_CHG.cube.ref OUT.autotest/SPIN1_CHG.cube 8
    echo "SPIN1_CHG.cube_pass $?" >>$1
    python3 $COMPARE_SCRIPT SPIN1_TAU.cube.ref OUT.autotest/SPIN1_TAU.cube 8
    echo "SPIN1_TAU.cube_pass $?" >>$1
fi

#---------------------------------------
# local density of states
# echo $has_ldos
#---------------------------------------
if ! test -z "$has_ldos"  && [  $has_ldos == 1 ]; then
    stm_bias=$(get_input_key_value "stm_bias" "OUT.autotest/INPUT")
    python3 $COMPARE_SCRIPT LDOS.cube.ref OUT.autotest/LDOS_"$stm_bias"eV.cube 8
    echo "LDOS.cube_pass $?" >> $1
fi

#---------------------------------------
# wave functions in real space 
# echo "$has_wfc_r" ## test out_wfc_r > 0
#---------------------------------------
if ! test -z "$has_wfc_r"  && [ $has_wfc_r == 1 ]; then
	if [[ ! -f OUT.autotest/running_scf.log ]];then
		echo "Can't find file OUT.autotest/running_scf.log"
		exit 1
	fi
	nband=$(grep NBANDS OUT.autotest/running_scf.log|awk '{print $3}')
    allgrid=$(grep "fft grid for wave functions" OUT.autotest/running_scf.log | awk -F "[=,\\\[\\\]]" '{print $3*$4*$5}')
	for((band=0;band<$nband;band++));do
		if [[ -f "OUT.autotest/wfc_realspace/wfc_realspace_0_$band" ]];then
			variance_wfc_r=`sed -n "13,$"p OUT.autotest/wfc_realspace/wfc_realspace_0_$band | \
						awk -v all=$allgrid 'BEGIN {sumall=0} {for(i=1;i<=NF;i++) {sumall+=($i-1)*($i-1)}}\
						END {printf"%.5f",(sumall/all)}'`
			echo "variance_wfc_r_0_$band $variance_wfc_r" >>$1
		else
			echo "Can't find file OUT.autotest/wfc_realspace/wfc_realspace_0_$band"
			exit 1
		fi
	done
fi	

#--------------------------------------------
# wave functions in plane wave basis 
# echo "$has_wfc_pw" ## test out_wfc_pw > 0
#--------------------------------------------
if ! test -z "$has_wfc_pw"  && [ $has_wfc_pw == 1 ]; then
	if [[ ! -f OUT.autotest/wfs1k1_pw.txt ]];then
		echo "Can't find file OUT.autotest/wfs1k1_pw.txt"
		exit 1
	fi
	awk 'BEGIN {max=0;read=0;band=1}
	{
		if(read==0 && $2 == "Band" && $3 == band){read=1}
		else if(read==1 && $2 == "Band" && $3 == band)
			{printf"Max_wfc_%d %.4f\n",band,max;read =0;band+=1;max=0}
		else if(read==1)
			{
				for(i=1;i<=NF;i++) 
				{
					if(sqrt($i*$i)>max) {max=sqrt($i*$i)}
				}
			} 
	}' OUT.autotest/wfs1k1_pw.txt >> $1
fi


#--------------------------------------------
# wave functions in LCAO basis
# echo "$has_lowf" # test out_wfc_lcao > 0
#--------------------------------------------
if ! test -z "$has_lowf"  && [ $has_lowf == 1 ]; then
	if ! test -z "$gamma_only"  && [ $gamma_only == 1 ]; then
		wfc_cal=OUT.autotest/wfs1_nao.txt
		wfc_ref=wfs1_nao.txt.ref
	else  # multi-k point case
		if ! test -z "$out_app_flag"  && [ $out_app_flag == 0 ]; then
			wfc_name=wfs1k1g3_nao
		else
			wfc_name=wfs1k2_nao
		fi
		awk 'BEGIN {flag=999}
    	{
        	if($2 == "(band)") {flag=2;print $0}
        	else if(flag>0) {flag-=1;print $0}
        	else if(flag==0) 
        	{
            	for(i=1;i<=NF/2;i++)
            	{printf "%.10e ",sqrt( $(2*i)*$(2*i)+$(2*i-1)*$(2*i-1) )};
            	printf "\n"
        	}	
        	else {print $0}
    	}' OUT.autotest/"$wfc_name".txt > OUT.autotest/"$wfc_name"_mod.txt
		wfc_cal=OUT.autotest/"$wfc_name"_mod.txt
		wfc_ref="$wfc_name"_mod.txt.ref
	fi

	python3 $COMPARE_SCRIPT $wfc_cal $wfc_ref 8 -abs 1
	echo "Compare_wfc_lcao_pass $?" >>$1
fi

#--------------------------------------------
# density matrix information 
#--------------------------------------------
if ! test -z "$out_dm"  && [ $out_dm == 1 ]; then
      dmfile=OUT.autotest/dms1_nao.txt
	  dmref=dms1_nao.txt.ref
      if test -z "$dmfile"; then
              echo "Can't find DM files"
              exit 1
      else
			python3 $COMPARE_SCRIPT $dmref $dmfile 5
            echo "DM_different $?" >>$1
      fi
fi

#--------------------------------------------
# mulliken charge
#--------------------------------------------
if ! test -z "$out_mul"  && [ $out_mul == 1 ]; then
    python3 $COMPARE_SCRIPT mulliken.txt.ref OUT.autotest/mulliken.txt 3
	echo "Compare_mulliken_pass $?" >>$1
fi

#--------------------------------------------
# Process .cube files for:
# 1. get_wf/get_pchg calculation tag (LCAO)
# 2. out_wfc_norm/out_wfc_re_im/out_pchg (PW)
#--------------------------------------------
need_process_cube=false
# Check if this is a LCAO calculation with get_wf/get_pchg
if [ $calculation == "get_wf" ] || [ $calculation == "get_pchg" ]; then
    need_process_cube=true
fi
# Check if this is a PW calculation with out_wfc_norm/out_wfc_re_im
out_wfc_norm=$(get_input_key_value "out_wfc_norm" "INPUT")
out_wfc_re_im=$(get_input_key_value "out_wfc_re_im" "INPUT")
out_pchg=$(get_input_key_value "out_pchg" "INPUT")
if [ -n "$out_wfc_norm" ] || [ -n "$out_wfc_re_im" ] || [ -n "$out_pchg" ]; then
    need_process_cube=true
fi
# Process .cube files if needed
if [ "$need_process_cube" = true ]; then
    cubefiles=$(ls OUT.autotest/ | grep -E '.cube$')
    
    if [ -z "$cubefiles" ]; then
        echo "Error: No .cube files found in OUT.autotest/"
        exit 1
    else
        for cube in $cubefiles; do
            total_chg=$($SUM_CUBE_EXE OUT.autotest/$cube)
            echo "$cube $total_chg" >> $1
        done
    fi
fi

#--------------------------------------------
# implicit solvation model
#--------------------------------------------
if ! test -z "$imp_sol" && [ $imp_sol == 1 ]; then
	esol_el=`grep E_sol_el $running_path | awk '{print $3}'`
	esol_cav=`grep E_sol_cav $running_path | awk '{print $3}'`
	echo "esolelref $esol_el" >>$1
	echo "esolcavref $esol_cav" >>$1
fi

#--------------------------------------------
# random phase approximation
#--------------------------------------------
if ! test -z "$run_rpa" && [ $run_rpa == 1 ]; then
	Etot_without_rpa=`grep Etot_without_rpa log.txt | awk 'BEGIN{FS=":"} {print $2}' `
	echo "Etot_without_rpa $Etot_without_rpa" >> $1
	onref=refcoulomb_mat_0.txt
	oncal=coulomb_mat_0.txt
	python3 $COMPARE_SCRIPT $onref $oncal 8
fi

#--------------------------------------------
# deepks
#--------------------------------------------
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
bash ${script_dir}/catch_deepks_properties.sh $1

#--------------------------------------------
# check symmetry 
#--------------------------------------------
if ! test -z "$symmetry" && [ $symmetry == 1 ]; then
	pointgroup=`grep 'POINT GROUP' $running_path | tail -n 2 | head -n 1 | awk '{print $4}'`
	spacegroup=`grep 'SPACE GROUP' $running_path | tail -n 1 | awk '{print $7}'`
	nksibz=`grep ' nkstot_ibz ' $running_path | awk '{print $3}'`
	echo "pointgroupref $pointgroup" >>$1
	echo "spacegroupref $spacegroup" >>$1
	echo "nksibzref $nksibz" >>$1
fi

#--------------------------------------------
# check currents in rt-TDDFT 
#--------------------------------------------
if ! test -z "$out_current" && [ $out_current ]; then
	current1ref=refcurrent_total.dat
	current1cal=OUT.autotest/current_total.dat
	python3 $COMPARE_SCRIPT $current1ref $current1cal 10
	echo "CompareCurrent_pass $?" >>$1
fi

#--------------------------------------------
# Linear response function 
#--------------------------------------------
if [ $is_lr == 1 ]; then
	lrns=$(get_input_key_value "lr_nstates" "INPUT")
	lrns1=`echo "$lrns + 1" |bc`
	grep -A$lrns1 "Excitation Energy" $running_path | awk 'NR > 2 && $2 ~ /^[0-9]+\.[0-9]+$/ {print $2}' > lr_eig.txt
	lreig_tot=`sum_file lr_eig.txt`
	echo "totexcitationenergyref $lreig_tot" >>$1
fi

#--------------------------------------------
# Check RDMFT method 
#--------------------------------------------
if ! test -z "$rdmft" && [[ $rdmft == 1 ]]; then
	echo "" >>$1
	echo "The following energy units are in Rydberg:" >>$1

	E_TV_RDMFT=$(grep "E_TV_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_TV_RDMFT_ref $E_TV_RDMFT" >>$1

	E_hartree_RDMFT=$(grep "E_hartree_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_hartree_RDMFT_ref $E_hartree_RDMFT" >>$1

	Exc_cwp22_RDMFT=$(grep "Exc_cwp22_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "Exc_cwp22_RDMFT_ref $Exc_cwp22_RDMFT" >>$1

	E_Ewald=$(grep "E_Ewald" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_Ewald_ref $E_Ewald" >>$1

	E_entropy=$(grep "E_entropy(-TS)" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_entropy_ref $E_entropy" >>$1

	E_descf=$(grep "E_descf" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_descf_ref $E_descf" >>$1

	Etotal_RDMFT=$(grep "Etotal_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "Etotal_RDMFT_ref $Etotal_RDMFT" >>$1

	Exc_ksdft=$(grep "Exc_ksdft" "$running_path" | tail -1 | awk '{print $2}')
	echo "Exc_ksdft_ref $Exc_ksdft" >>$1

	E_exx_ksdft=$(grep "E_exx_ksdft" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_exx_ksdft_ref $E_exx_ksdft" >>$1

	echo "" >>$1
fi

#--------------------------------------------
# Check time information 
#--------------------------------------------
#echo $total_band
ttot=`grep $word_total_time $running_path | awk '{print $3}'`
echo "totaltimeref $ttot" >>$1
