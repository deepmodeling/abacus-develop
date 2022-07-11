#!/bin/bash

sum_file(){
	inc=1
	if ! test -z $2; then
		inc=$2
	fi

    awk -v each=$inc -F '[ (),]' 'BEGIN {sum=0.0;line=1} 
        {
			if(NR==line) 
		 	{
				for(i=1;i<=NF;i++) {sum+=sqrt($i*$i)};
				line+=each
			}
		}
        END {printf "%.6f", sum}' $1
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
has_hs=`grep -En '(^|[[:space:]])out_mat_hs($|[[:space:]])' INPUT | awk '{print $2}'`
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
#echo $running_path
base=`grep -En '(^|[[:space:]])basis_type($|[[:space:]])' INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
if [ $base == "pw" ]; then word="plane_wave_line" 
else
word="lcao_line"
fi
#echo $word
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
if ! test -z "$has_force" && [ $has_force -eq 1 ]; then
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
if ! test -z "$has_stress" && [  $has_stress -eq 1 ]; then
	#grep -A6 "TOTAL-STRESS" $running_path|sed '1,4d'|sed '4,8d' >stress.txt
    grep -A6 "TOTAL-STRESS" $running_path| awk 'NF==3' | tail -3> stress.txt
	total_stress=`sum_file stress.txt`
	rm stress.txt
	echo "totalstressref $total_stress" >>$1
fi


#echo $total_stress
#if ! test -z "$has_charge" && [  $has_charge -eq 1 ]; then
#	total_charge=`sum_file OUT.autotest/SPIN1_CHG`
#	echo "totalchargeref $total_charge" >>$1
#fi


#echo $total_charge
if ! test -z "$has_dos"  && [  $has_dos -eq 1 ]; then
	total_dos=`cat OUT.autotest/DOS1_smearing.dat | awk 'END {print}' | awk '{print $3}'`
	echo "totaldosref $total_dos" >> $1
fi
#	smearing_dos=`sum_file OUT.autotest/DOS1_smearing.dat`
#	echo "totaldossmearing $smearing_dos" >> $1


#echo total_dos
#echo $has_band
if ! test -z "$has_band"  && [  $has_band -eq 1 ]; then
	total_band=`sum_file OUT.autotest/BANDS_1.dat`
	echo "totalbandref $total_band" >>$1
fi
#echo $has_hs
if ! test -z "$has_hs"  && [  $has_hs -eq 1 ]; then
        total_h=`sum_file OUT.autotest/data-0-H`
        echo "totalHmatrix $total_h" >>$1
	total_s=`sum_file OUT.autotest/data-0-S`
	echo "totalSmatrix $total_s" >>$1
fi

# echo "$has_wfc_r" ## test out_wfc_r > 0
if ! test -z "$has_wfc_r"  && [ $has_wfc_r -eq 1 ]; then
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
if ! test -z "$has_wfc_pw"  && [ $has_wfc_pw -eq 1 ]; then
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
if ! test -z "$has_lowf"  && [ $has_lowf -eq 1 ]; then
	if ! test -z "$gamma_only"  && [ $gamma_only -eq 1 ]; then
		lowfiles=`ls OUT.autotest/ | grep LOWF_GAMMA`
	else
		lowfiles=`ls OUT.autotest/ | grep LOWF_K`
	fi
	if test -z "$lowfiles"; then
		echo "Can't find LOWF files"
		exit 1
	else
		for lowf in $lowfiles;
		do
			if ! test -f OUT.autotest/$lowf; then
				echo "Irregular LOWF file found"
				exit 1
			else
				sed -i "1,$ s/[a-d]//g" OUT.autotest/$lowf
				sed -i "1,$ s/[f-z]//g" OUT.autotest/$lowf
				sed -i "1,$ s/[A-D]//g" OUT.autotest/$lowf
				sed -i "1,$ s/[F-Z]//g" OUT.autotest/$lowf
				sed -i "1,$ s/)//g" OUT.autotest/$lowf
				sed -i "1,$ s/(//g" OUT.autotest/$lowf
				total_lowf=`sum_file OUT.autotest/$lowf`
				echo "$lowf $total_lowf" >>$1
			fi
		done
	fi
fi

if ! test -z "$out_dm"  && [ $out_dm -eq 1 ]; then
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
                              sed -i "1,$ s/[a-d]//g" OUT.autotest/$dm
                              sed -i "1,$ s/[f-z]//g" OUT.autotest/$dm
                              sed -i "1,$ s/[A-D]//g" OUT.autotest/$dm
                              sed -i "1,$ s/[F-Z]//g" OUT.autotest/$dm
                              sed -i "1,$ s/)//g" OUT.autotest/$dm
                              sed -i "1,$ s/(//g" OUT.autotest/$dm
                              total_dm=`sum_file OUT.autotest/$dm`
                              echo "$dm $total_dm" >>$1
                      fi
              done
      fi
fi

if ! test -z "$out_mul"  && [ $out_mul -eq 1 ]; then
      mulfile=`ls OUT.autotest/ | grep mulliken`
      if test -z "$mulfile"; then
              echo "Can't find Mulliken files"
              exit 1
      else
              for mul in $mulfile;
              do
                      if ! test -f OUT.autotest/$mul; then
                              echo "Irregular Mulliken file found"
                              exit 1
                      else
                              sed -i "1,$ s/[a-d]//g" OUT.autotest/$mul
                              sed -i "1,$ s/[f-z]//g" OUT.autotest/$mul
                              sed -i "1,$ s/[A-D]//g" OUT.autotest/$mul
                              sed -i "1,$ s/[F-Z]//g" OUT.autotest/$mul
                              sed -i "1,$ s/+//g" OUT.autotest/$mul
                              sed -i "1,$ s/)//g" OUT.autotest/$mul
                              sed -i "1,$ s/(//g" OUT.autotest/$mul
                              total_mul=`sum_file OUT.autotest/$mul`
                              echo "$mul $total_mul" >>$1
			fi
		done
	fi
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

if ! test -z "$imp_sol" && [ $imp_sol -eq 1 ]; then
	esol_el=`grep E_sol_el $running_path | awk '{print $3}'`
	esol_cav=`grep E_sol_cav $running_path | awk '{print $3}'`
	echo "esolelref $esol_el" >>$1
	echo "esolcavref $esol_cav" >>$1
fi

#echo $total_band
ttot=`grep $word $running_path | awk '{print $3}'`
echo "totaltimeref $ttot" >>$1

if ! test -z "$deepks_out_labels" && [ $deepks_out_labels -eq 1 ]; then
	sed '/n_des/d' descriptor.dat > des_tmp.txt
	total_des=`sum_file des_tmp.txt 5`
	rm des_tmp.txt
	echo "totaldes $total_des" >>$1
fi

if ! test -z "$deepks_bandgap" && [ $deepks_bandgap -eq 1 ]; then
	odelta=`python3 get_odelta.py`
	echo "odelta $odelta" >>$1
	oprec=`python3 get_oprec.py`
	echo "oprec $oprec" >> $1
fi
