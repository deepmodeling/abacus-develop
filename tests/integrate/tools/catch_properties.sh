#!/bin/bash

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
#answer=`sum_file test.txt`
#echo $answer
#exit 0

get_input_key_value(){
	key=$1
	inputf=$2
	value=$(awk -v key=$key '{if($1==key) a=$2} END {print a}' $inputf)
	echo $value
}

file=$1
#echo $1
calculation=`grep calculation INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
running_path=`echo "OUT.autotest/running_$calculation"".log"`
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
has_eband_separate=$(get_input_key_value "out_eband_terms" "INPUT")
has_r=$(get_input_key_value "out_mat_r" "INPUT")
deepks_out_labels=$(get_input_key_value "deepks_out_labels" "INPUT")
deepks_scf=$(get_input_key_value "deepks_scf" "INPUT")
deepks_bandgap=$(get_input_key_value "deepks_bandgap" "INPUT")
deepks_v_delta=$(get_input_key_value "deepks_v_delta" "INPUT")
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
esolver_type=$(get_input_key_value "esolver_type" "INPUT")
#echo $running_path
base=$(get_input_key_value "basis_type" "INPUT")
word="driver_line"
symmetry=$(get_input_key_value "symmetry" "INPUT")
out_current=$(get_input_key_value "out_current" "INPUT")
test -e $1 && rm $1
#--------------------------------------------
# if NOT non-self-consistent calculations or linear response
#--------------------------------------------
is_lr=0
if [ ! -z $esolver_type ] && ([ $esolver_type == "lr" ] || [ $esolver_type == "ks-lr" ]); then
	is_lr=1
fi

if [ $calculation != "nscf" ] && [ $calculation != "get_wf" ]\
&& [ $calculation != "get_pchg" ] && [ $calculation != "get_S" ]\
&& [ $is_lr == 0 ]; then
	#etot=`grep ETOT_ $running_path | awk '{print $2}'` 
	etot=$(grep "ETOT_" "$running_path" | tail -1 | awk '{print $2}')
	etotperatom=`awk 'BEGIN {x='$etot';y='$natom';printf "%.10f\n",x/y}'`
	echo "etotref $etot" >>$1
	echo "etotperatomref $etotperatom" >>$1
fi


#echo $etot
#echo "hasforce:"$has_force
if ! test -z "$has_force" && [ $has_force == 1 ]; then
	nn3=`echo "$natom + 1" |bc`
	#nn1=`echo "$natom + 1" |bc`
	#nn5=`echo "$natom + 6" |bc`
	#grep -A$nn3 "TOTAL-FORCE" $running_path|sed '1,5d'|sed ''$nn1','$nn5'd'|awk '{printf $2"\t"$3"\t"$4"\n"}' > force.txt
    grep -A$nn3 "TOTAL-FORCE" $running_path |awk 'NF==4{print $2,$3,$4}' | tail -$natom > force.txt  #check the last step result
	total_force=`sum_file force.txt`
	rm force.txt
	echo "totalforceref $total_force" >>$1
fi

#echo $total_force
#echo "has_stress:"$has_stress
if ! test -z "$has_stress" && [  $has_stress == 1 ]; then
	#grep -A6 "TOTAL-STRESS" $running_path|sed '1,4d'|sed '4,8d' >stress.txt
    grep -A4 "TOTAL-STRESS" $running_path| awk 'NF==3' | tail -3> stress.txt
	total_stress=`sum_file stress.txt`
	rm stress.txt
	echo "totalstressref $total_stress" >>$1
fi


#echo $total_stress
#if ! test -z "$has_charge" && [  $has_charge == 1 ]; then
#	total_charge=`sum_file OUT.autotest/SPIN1_CHG`
#	echo "totalchargeref $total_charge" >>$1
#fi


#echo $total_charge
if ! test -z "$has_dos"  && [  $has_dos == 1 ]; then
	total_dos=`cat OUT.autotest/DOS1_smearing.dat | awk 'END {print}' | awk '{print $3}'`
	echo "totaldosref $total_dos" >> $1
fi
#	smearing_dos=`sum_file OUT.autotest/DOS1_smearing.dat`
#	echo "totaldossmearing $smearing_dos" >> $1

#echo Onsager coefficiency
if ! test -z "$has_cond"  && [  $has_cond == 1 ]; then
	onref=refOnsager.txt
	oncal=Onsager.txt
	python3 ../tools/CompareFile.py $onref $oncal 3 -com_type 0
    echo "CompareH_Failed $?" >>$1
	rm -f je-je.txt Chebycoef
fi

#echo $out_dm1
if ! test -z "$out_dm1"  && [  $out_dm1 == 1 ]; then
	dm1ref=refdata-DMR-sparse_SPIN0.csr
	dm1cal=OUT.autotest/data-DMR-sparse_SPIN0.csr
	python3 ../tools/CompareFile.py $dm1ref $dm1cal 8
	echo "CompareDM1_pass $?" >>$1
fi

#echo $out_pot1
if ! test -z "$out_pot"  && [  $out_pot == 1 ]; then
	pot1ref=refSPIN1_POT.cube
	pot1cal=OUT.autotest/SPIN1_POT.cube
	python3 ../tools/CompareFile.py $pot1ref $pot1cal 3
	echo "ComparePot1_pass $?" >>$1
fi

#echo $out_pot2
if ! test -z "$out_pot"  && [  $out_pot == 2 ]; then
	pot1ref=refElecStaticPot.cube
	pot1cal=OUT.autotest/ElecStaticPot.cube
	python3 ../tools/CompareFile.py $pot1ref $pot1cal 8
	echo "ComparePot_pass $?" >>$1
fi

#echo $out_elf
if ! test -z "$out_elf"  && [  $out_elf == 1 ]; then
	elf1ref=refELF.cube
	elf1cal=OUT.autotest/ELF.cube
	python3 ../tools/CompareFile.py $elf1ref $elf1cal 3
	echo "ComparePot1_pass $?" >>$1
fi

#echo $get_s
if ! test -z "$get_s"  && [  $get_s == "get_S" ]; then
	sref=refSR.csr
	scal=OUT.autotest/SR.csr
	python3 ../tools/CompareFile.py $sref $scal 8
	echo "CompareS_pass $?" >>$1
fi

#echo $out_pband
if ! test -z "$out_pband"  && [  $out_pband == 1 ]; then
	#pbandref=refPBANDS_1
	#pbandcal=OUT.autotest/PBANDS_1
	#python3 ../tools/CompareFile.py $pbandref $pbandcal 8
	#echo "CompareProjBand_pass $?" >>$1
	orbref=refOrbital
	orbcal=OUT.autotest/Orbital
	python3 ../tools/CompareFile.py $orbref $orbcal 8
	echo "CompareOrb_pass $?" >>$1
fi

#echo $toW90
if ! test -z "$toW90"  && [  $toW90 == 1 ]; then
	amnref=diamond.amn
	amncal=OUT.autotest/diamond.amn
	mmnref=diamond.mmn
	mmncal=OUT.autotest/diamond.mmn
	eigref=diamond.eig
	eigcal=OUT.autotest/diamond.eig
	sed -i '1d' $amncal
	sed -i '1d' $mmncal
	python3 ../tools/CompareFile.py $amnref $amncal 1 -abs 8
	echo "CompareAMN_pass $?" >>$1
	python3 ../tools/CompareFile.py $mmnref $mmncal 1 -abs 8
	echo "CompareMMN_pass $?" >>$1
	python3 ../tools/CompareFile.py $eigref $eigcal 8
	echo "CompareEIG_pass $?" >>$1
fi

#echo total_dos
#echo $has_band
if ! test -z "$has_band"  && [  $has_band == 1 ]; then
	bandref=refBANDS_1.dat
	bandcal=OUT.autotest/BANDS_1.dat
	python3 ../tools/CompareFile.py $bandref $bandcal 8
	echo "CompareBand_pass $?" >>$1
fi
#echo $has_hs
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

        python3 ../tools/CompareFile.py $href $hcal 6
    echo "CompareH_pass $?" >>$1
    python3 ../tools/CompareFile.py $sref $scal 8
    echo "CompareS_pass $?" >>$1
fi

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
	python3 ../tools/CompareFile.py $xcref $xccal 4
	echo "CompareVXC_pass $?" >>$1
	python3 ../tools/CompareFile.py $oeref $oecal 5
    echo "CompareOrbXC_pass $?" >>$1
fi

if ! test -z "$has_eband_separate"  && [  $has_eband_separate == 1 ]; then
	ekref=kinetic_out.ref
	ekcal=OUT.autotest/kinetic_out.dat
	python3 ../tools/CompareFile.py $ekref $ekcal 4
	echo "CompareOrbKinetic_pass $?" >>$1
	vlref=vpp_local_out.ref
	vlcal=OUT.autotest/vpp_local_out.dat
	python3 ../tools/CompareFile.py $vlref $vlcal 4
	echo "CompareOrbVL_pass $?" >>$1
	vnlref=vpp_nonlocal_out.ref
	vnlcal=OUT.autotest/vpp_nonlocal_out.dat
	python3 ../tools/CompareFile.py $vnlref $vnlcal 4
	echo "CompareOrbVNL_pass $?" >>$1
	vhref=vhartree_out.ref
	vhcal=OUT.autotest/vhartree_out.dat
	python3 ../tools/CompareFile.py $vhref $vhcal 4
	echo "CompareOrbVHartree_pass $?" >>$1
fi

#echo $has_hs2
if ! test -z "$has_hs2"  && [  $has_hs2 == 1 ]; then
    #python3 ../tools/CompareFile.py data-HR-sparse_SPIN0.csr.ref OUT.autotest/data-HR-sparse_SPIN0.csr 8
    #echo "CompareHR_pass $?" >>$1
    python3 ../tools/CompareFile.py data-SR-sparse_SPIN0.csr.ref OUT.autotest/data-SR-sparse_SPIN0.csr 8
    echo "CompareSR_pass $?" >>$1
fi

#echo $has_mat_r
if ! test -z "$has_mat_r"  && [  $has_mat_r == 1 ]; then
    python3 ../tools/CompareFile.py data-rR-sparse.csr.ref OUT.autotest/data-rR-sparse.csr 8
    echo "ComparerR_pass $?" >>$1
fi

#echo $has_mat_t
if ! test -z "$has_mat_t"  && [  $has_mat_t == 1 ]; then
    python3 ../tools/CompareFile.py data-TR-sparse_SPIN0.csr.ref OUT.autotest/data-TR-sparse_SPIN0.csr 8
    echo "ComparerTR_pass $?" >>$1
fi

#echo $has_mat_dh
if ! test -z "$has_mat_dh"  && [  $has_mat_dh == 1 ]; then
    python3 ../tools/CompareFile.py data-dHRx-sparse_SPIN0.csr.ref OUT.autotest/data-dHRx-sparse_SPIN0.csr 8
    echo "ComparerdHRx_pass $?" >>$1
    python3 ../tools/CompareFile.py data-dHRy-sparse_SPIN0.csr.ref OUT.autotest/data-dHRy-sparse_SPIN0.csr 8
    echo "ComparerdHRy_pass $?" >>$1
    python3 ../tools/CompareFile.py data-dHRz-sparse_SPIN0.csr.ref OUT.autotest/data-dHRz-sparse_SPIN0.csr 8
    echo "ComparerdHRz_pass $?" >>$1
fi

#echo $has_scan
if ! test -z "$has_scan"  && [  $has_scan == "scan" ] && \
       ! test -z "$out_chg" && [ $out_chg == 1 ]; then
    python3 ../tools/CompareFile.py SPIN1_CHG.cube.ref OUT.autotest/SPIN1_CHG.cube 8
    echo "SPIN1_CHG.cube_pass $?" >>$1
    python3 ../tools/CompareFile.py SPIN1_TAU.cube.ref OUT.autotest/SPIN1_TAU.cube 8
    echo "SPIN1_TAU.cube_pass $?" >>$1
fi

# echo "$has_wfc_r" ## test out_wfc_r > 0
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

# echo "$has_wfc_pw" ## test out_wfc_pw > 0
if ! test -z "$has_wfc_pw"  && [ $has_wfc_pw == 1 ]; then
	if [[ ! -f OUT.autotest/WAVEFUNC1.txt ]];then
		echo "Can't find file OUT.autotest/WAVEFUNC1.txt"
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
	}' OUT.autotest/WAVEFUNC1.txt >> $1
fi

# echo "$has_lowf" ## test out_wfc_lcao > 0
if ! test -z "$has_lowf"  && [ $has_lowf == 1 ]; then
	if ! test -z "$gamma_only"  && [ $gamma_only == 1 ]; then
		wfc_cal=OUT.autotest/WFC_NAO_GAMMA1.txt
		wfc_ref=WFC_NAO_GAMMA1.txt.ref	
	else
		if ! test -z "$out_app_flag"  && [ $out_app_flag == 0 ]; then
			wfc_name=WFC_NAO_K1_ION3
		else
			wfc_name=WFC_NAO_K2
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

	python3 ../tools/CompareFile.py $wfc_cal $wfc_ref 8 -abs 1
	echo "Compare_wfc_lcao_pass $?" >>$1
fi

if ! test -z "$out_dm"  && [ $out_dm == 1 ]; then
      dmfile=OUT.autotest/SPIN1_DM
	  dmref=SPIN1_DM.ref
      if test -z "$dmfile"; then
              echo "Can't find DM files"
              exit 1
      else
			python3 ../tools/CompareFile.py $dmref $dmfile 5
            echo "DM_different $?" >>$1
      fi
fi

if ! test -z "$out_mul"  && [ $out_mul == 1 ]; then
    python3 ../tools/CompareFile.py mulliken.txt.ref OUT.autotest/mulliken.txt 6
	echo "Compare_mulliken_pass $?" >>$1
fi

if [ $calculation == "get_wf" ]; then
	nfile=0
	# envfiles=`ls OUT.autotest/ | grep ENV$`
	# if test -z "$envfiles"; then
	# 	echo "Can't find ENV(-elope) files"
	# 	exit 1
	# else
	# 	for env in $envfiles;
	# 	do
	# 		nelec=`../tools/sum_ENV_H2 OUT.autotest/$env`
	# 		nfile=$(($nfile+1))
	# 		echo "nelec$nfile $nelec" >>$1
	# 	done
	# fi
	cubefiles=`ls OUT.autotest/ | grep -E '.cube$'`
	if test -z "$cubefiles"; then
		echo "Can't find BAND_CHG files"
		exit 1
	else
		for cube in $cubefiles;
		do
			total_chg=`../tools/sum_ENV_H2_cube OUT.autotest/$cube`
			echo "$cube $total_chg" >>$1
		done
	fi
fi

if [ $calculation == "get_pchg" ]; then
	nfile=0
	# chgfiles=`ls OUT.autotest/ | grep -E '_CHG$'`
	# if test -z "$chgfiles"; then
	# 	echo "Can't find BAND_CHG files"
	# 	exit 1
	# else
	# 	for chg in $chgfiles;
	# 	do
	# 		total_chg=`../tools/sum_BAND_CHG_H2 OUT.autotest/$chg`
	# 		nfile=$(($nfile+1))
	# 		echo "nelec$nfile $total_chg" >>$1
	# 	done
	# fi
	cubefiles=`ls OUT.autotest/ | grep -E '.cube$'`
	if test -z "$cubefiles"; then
		echo "Can't find BAND_CHG files"
		exit 1
	else
		for cube in $cubefiles;
		do
			total_chg=`../tools/sum_BAND_CHG_H2_cube OUT.autotest/$cube`
			echo "$cube $total_chg" >>$1
		done
	fi
fi

if ! test -z "$imp_sol" && [ $imp_sol == 1 ]; then
	esol_el=`grep E_sol_el $running_path | awk '{print $3}'`
	esol_cav=`grep E_sol_cav $running_path | awk '{print $3}'`
	echo "esolelref $esol_el" >>$1
	echo "esolcavref $esol_cav" >>$1
fi

if ! test -z "$run_rpa" && [ $run_rpa == 1 ]; then
	Etot_without_rpa=`grep Etot_without_rpa log.txt | awk 'BEGIN{FS=":"} {print $2}' `
	echo "Etot_without_rpa $Etot_without_rpa" >> $1
	onref=refcoulomb_mat_0.txt
	oncal=coulomb_mat_0.txt
	python3 ../tools/CompareFile.py $onref $oncal 8
fi

if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 1 ]; then
	sed '/n_des/d' OUT.autotest/deepks_desc.dat > des_tmp.txt
	total_des=`sum_file des_tmp.txt 5`
	rm des_tmp.txt
	echo "totaldes $total_des" >>$1
	if ! test -z "$deepks_scf" && [ $deepks_scf == 1 ]; then
		deepks_e_dm=`python3 get_dm_eig.py`
	    echo "deepks_e_dm $deepks_e_dm" >>$1
	fi
	if ! test -z "$has_force" && [ $has_force == 1 ]; then
	    deepks_f_label=`python3 get_grad_vx.py`
		echo "deepks_f_label $deepks_f_label" >>$1
	fi
	if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
	    deepks_s_label=`python3 get_grad_vepsl.py`
		echo "deepks_s_label $deepks_s_label" >>$1
	fi
fi

if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ]; then
	odelta=`python3 get_odelta.py`
	echo "odelta $odelta" >>$1
	oprec=`python3 get_oprec.py`
	echo "oprec $oprec" >> $1
fi

if ! test -z "$deepks_v_delta" && [ $deepks_v_delta == 1 ]; then
	totalh=`python3 get_sum_numpy.py OUT.autotest/deepks_htot.npy `
	echo "totalh $totalh" >>$1
	totalvdelta=`python3 get_v_delta.py`
	echo "totalvdelta $totalvdelta" >>$1
	totalvdp=`python3 get_sum_numpy.py OUT.autotest/deepks_vdpre.npy `
	echo "totalvdp $totalvdp" >> $1
fi

if ! test -z "$deepks_v_delta" && [ $deepks_v_delta == 2 ]; then
	totalh=`python3 get_sum_numpy.py OUT.autotest/deepks_htot.npy `
	echo "totalh $totalh" >>$1
	totalvdelta=`python3 get_v_delta.py`
	echo "totalvdelta $totalvdelta" >>$1
	total_psialpha=`python3 get_sum_numpy.py OUT.autotest/deepks_psialpha.npy `
	echo "total_psialpha $total_psialpha" >> $1
	total_gevdm=`python3 get_sum_numpy.py OUT.autotest/deepks_gevdm.npy `
	echo "total_gevdm $total_gevdm" >> $1
fi

if ! test -z "$symmetry" && [ $symmetry == 1 ]; then
	pointgroup=`grep 'POINT GROUP' $running_path | tail -n 2 | head -n 1 | awk '{print $4}'`
	spacegroup=`grep 'SPACE GROUP' $running_path | tail -n 1 | awk '{print $7}'`
	nksibz=`grep ' nkstot_ibz ' $running_path | awk '{print $3}'`
	echo "pointgroupref $pointgroup" >>$1
	echo "spacegroupref $spacegroup" >>$1
	echo "nksibzref $nksibz" >>$1
fi

if ! test -z "$out_current" && [ $out_current ]; then
	current1ref=refcurrent_total.dat
	current1cal=OUT.autotest/current_total.dat
	python3 ../tools/CompareFile.py $current1ref $current1cal 10
	echo "CompareCurrent_pass $?" >>$1
fi

if [ $is_lr == 1 ]; then
	lrns=$(get_input_key_value "lr_nstates" "INPUT")
	lrns1=`echo "$lrns + 1" |bc`
	grep -A$lrns1 "Excitation Energy" $running_path | awk 'NR > 2 && $2 ~ /^[0-9]+\.[0-9]+$/ {print $2}' > lr_eig.txt
	lreig_tot=`sum_file lr_eig.txt`
	echo "totexcitationenergyref $lreig_tot" >>$1
fi

#echo $total_band
ttot=`grep $word $running_path | awk '{print $3}'`
echo "totaltimeref $ttot" >>$1
