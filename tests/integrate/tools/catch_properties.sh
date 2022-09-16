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

file=$1
#echo $1
calculation=`grep calculation INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
running_path=`echo "OUT.autotest/running_$calculation"".log"`
natom=`grep -En '(^|[[:space:]])TOTAL ATOM NUMBER($|[[:space:]])' $running_path | awk '{print $6}'`
has_force=`grep -En '(^|[[:space:]])cal_force($|[[:space:]])' INPUT | awk '{print $2}'`
has_stress=`grep -En '(^|[[:space:]])cal_stress($|[[:space:]])' INPUT | awk '{print $2}'`
has_dftu=`grep -En '(^|[[:space:]])dft_plus_u($|[[:space:]])' INPUT | awk '{print $2}'`
has_band=`grep -En '(^|[[:space:]])out_band($|[[:space:]])' INPUT | awk '{print $2}'`
has_dos=`grep -En '(^|[[:space:]])out_dos($|[[:space:]])' INPUT | awk '{print $2}'`
has_cond=`grep -En '(^|[[:space:]])cal_cond($|[[:space:]])' INPUT | awk '{print $2}'`
has_hs=`grep -En '(^|[[:space:]])out_mat_hs($|[[:space:]])' INPUT | awk '{print $2}'`
has_hs2=`grep -En '(^|[[:space:]])out_mat_hs2($|[[:space:]])' INPUT | awk '{print $2}'`
has_r=`grep -En '(^|[[:space:]])out_mat_r($|[[:space:]])' INPUT | awk '{print $2}'`
deepks_out_labels=`grep deepks_out_labels INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
deepks_bandgap=`grep deepks_bandgap INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
has_lowf=`grep out_wfc_lcao INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
has_wfc_r=`grep out_wfc_r INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
has_wfc_pw=`grep out_wfc_pw INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
out_dm=`grep out_dm INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
out_mul=`grep out_mul INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
gamma_only=`grep gamma_only INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
imp_sol=`grep imp_sol INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
comp_chg=`grep comp_chg INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
#echo $running_path
base=`grep -En '(^|[[:space:]])basis_type($|[[:space:]])' INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
word="driver_line"
test -e $1 && rm $1
#--------------------------------------------
# if NOT non-self-consistent calculations
#--------------------------------------------
if [ $calculation != "nscf" ] && [ $calculation != "ienvelope" ]\
&& [ $calculation != "istate" ]	; then
	etot=`grep ETOT_ $running_path | awk '{print $2}'`
	etotperatom=`awk 'BEGIN {x='$etot';y='$natom';printf "%.10f\n",x/y}'`
	echo "etotref $etot" >>$1
	echo "etotperatomref $etotperatom" >>$1
fi


#echo $etot
#echo "hasforce:"$has_force
if ! test -z "$has_force" && [ $has_force == 1 ]; then
	nn3=`echo "$natom + 4" |bc`
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
    grep -A6 "TOTAL-STRESS" $running_path| awk 'NF==3' | tail -3> stress.txt
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
	python3 ../tools/CompareFile.py $onref $oncal 2
    echo "CompareH_Failed $?" >>$1
	rm -f je-je.txt Chebycoef
fi

#echo total_dos
#echo $has_band
if ! test -z "$has_band"  && [  $has_band == 1 ]; then
	total_band=`sum_file OUT.autotest/BANDS_1.dat`
	echo "totalbandref $total_band" >>$1
fi
#echo $has_hs
if ! test -z "$has_hs"  && [  $has_hs == 1 ]; then
	if ! test -z "$gamma_only"  && [ $gamma_only == 1 ]; then
		href=data-0-H.ref
		hcal=OUT.autotest/data-0-H
		sref=data-0-H.ref
		scal=OUT.autotest/data-0-H
	else
		href=data-1-H.ref
		hcal=OUT.autotest/data-1-H
		sref=data-1-H.ref
		scal=OUT.autotest/data-1-H
	fi

	python3 ../tools/CompareFile.py $href $hcal 8
    echo "CompareH_pass $?" >>$1
    python3 ../tools/CompareFile.py $sref $scal 8
    echo "CompareS_pass $?" >>$1
fi

#echo $has_hs2
if ! test -z "$has_hs2"  && [  $has_hs2 == 1 ]; then
    python3 ../tools/CompareFile.py data-HR-sparse_SPIN0.csr.ref OUT.autotest/data-HR-sparse_SPIN0.csr 8
    echo "CompareHR_pass $?" >>$1
    python3 ../tools/CompareFile.py data-SR-sparse_SPIN0.csr.ref OUT.autotest/data-SR-sparse_SPIN0.csr 8
    echo "CompareSR_pass $?" >>$1
fi

# echo "$has_wfc_r" ## test out_wfc_r > 0
if ! test -z "$has_wfc_r"  && [ $has_wfc_r == 1 ]; then
	if [[ ! -f OUT.autotest/running_scf.log ]];then
		echo "Can't find file OUT.autotest/running_scf.log"
		exit 1
	fi
	nband=$(grep NBANDS OUT.autotest/running_scf.log|awk '{print $3}')
	allgrid=$(grep "fft grid for wave functions" OUT.autotest/running_scf.log|awk -F "[=,]" '{print $2*$3*$4}')
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
		wfc_cal=OUT.autotest/LOWF_GAMMA_S1.dat
		wfc_ref=LOWF_GAMMA_S1.dat.ref	
	else
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
    	}' OUT.autotest/LOWF_K_2.dat > OUT.autotest/LOWF_K_2_mod.dat
		wfc_cal=OUT.autotest/LOWF_K_2_mod.dat
		wfc_ref=LOWF_K_2_mod.dat.ref
	fi

	python3 ../tools/CompareFile.py $wfc_cal $wfc_ref 8 -abs 1
	echo "Compare_wfc_lcao_pass $?" >>$1
fi

if ! test -z "$out_dm"  && [ $out_dm == 1 ]; then
      dmfile=`ls OUT.autotest/ | grep "^SPIN1_DM"`
      if test -z "$dmfile"; then
              echo "Can't find DM files"
              exit 1
      else
              for dm in $dmfile;
              do
                      if ! test -f OUT.autotest/$dm; then
                              echo "Irregular DM file found"
                              exit 1
                      else
                            total_dm=$(awk 'BEGIN {sum=0.0;startline=999}
							{
								if(NR==7){startline=$1+14;}
								else if(NR>=startline) 
								{
									for(i=1;i<=NF;i++){sum+=sqrt($i*$i)}
								}
							}
							END {printf"%.6f",sum}' OUT.autotest/$dm)
                            echo "$dm $total_dm" >>$1
                      fi
              done
      fi
fi

if ! test -z "$out_mul"  && [ $out_mul == 1 ]; then
    python3 ../tools/CompareFile.py mulliken.txt.ref OUT.autotest/mulliken.txt 8
	echo "Compare_mulliken_pass $?" >>$1
fi

if [ $calculation == "ienvelope" ]; then
	nfile=0
	envfiles=`ls OUT.autotest/ | grep ENV`
	if test -z "$envfiles"; then
		echo "Can't find ENV(-elope) files"
		exit 1
	else
		for env in $envfiles;
		do
			nelec=`../tools/sum_ENV_H2 OUT.autotest/$env`
			nfile=$(($nfile+1))
			echo "nelec$nfile $nelec" >>$1	
		done
	fi
fi

if [ $calculation == "istate" ]; then
	chgfiles=`ls OUT.autotest/ | grep -E '_CHG$'`
	if test -z "$chgfiles"; then
		echo "Can't find BAND_CHG files"
		exit 1
	else
		for chg in $chgfiles;
		do
			total_chg=`../tools/sum_BAND_CHG_H2 OUT.autotest/$chg`
			echo "$chg $total_chg" >>$1
		done
	fi
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

if ! test -z "$comp_chg" && [ $comp_chg == 1 ]; then
	ecomp_self=`grep E_comp_self $running_path | awk '{print $3}'`
	ecomp_electron=`grep E_comp_electron $running_path | awk '{print $3}'`
	ecomp_nuclear=`grep E_comp_nuclear $running_path | awk '{print $3}'`
	ecomp_tot=`grep E_comp_tot $running_path | awk '{print $3}'`
	echo "ecompselfref $ecomp_self" >>$1
	echo "ecompelectronref $ecomp_electron" >>$1
	echo "ecompnuclearref $ecomp_nuclear" >>$1
	echo "ecomptotref $ecomp_tot" >>$1
fi

#echo $total_band
ttot=`grep $word $running_path | awk '{print $3}'`
echo "totaltimeref $ttot" >>$1

if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 1 ]; then
	sed '/n_des/d' descriptor.dat > des_tmp.txt
	total_des=`sum_file des_tmp.txt 5`
	rm des_tmp.txt
	echo "totaldes $total_des" >>$1
fi

if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ]; then
	odelta=`python3 get_odelta.py`
	echo "odelta $odelta" >>$1
	oprec=`python3 get_oprec.py`
	echo "oprec $oprec" >> $1
fi
