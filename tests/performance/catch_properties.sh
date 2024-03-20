#!/bin/bash

sum_file(){
	line=`grep -vc '^$' $1`
	sum=0.0
	for (( num=1 ; num<=$line ; num++ ));do
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
suffix=`awk '$1=="suffix"{print $2}' INPUT`
if [[ "$suffix" == "" ]];then
    suffix="ABACUS"
fi
calculation=`awk '$1=="calculation"{print $2}' INPUT`
if [[ "$suffix" == "" ]];then
    calculation="scf"
fi
running_path=`echo "OUT.${suffix}/running_$calculation"".log"`
natom=`grep -En '(^|[[:space:]])TOTAL ATOM NUMBER($|[[:space:]])' $running_path | awk '{print $6}'`
has_force=`grep -En '(^|[[:space:]])cal_force($|[[:space:]])' INPUT | awk '{print $2}'`
has_stress=`grep -En '(^|[[:space:]])cal_stress($|[[:space:]])' INPUT | awk '{print $2}'`
has_dftu=`grep -En '(^|[[:space:]])dft_plus_u($|[[:space:]])' INPUT | awk '{print $2}'`
has_band=`grep -En '(^|[[:space:]])out_band($|[[:space:]])' INPUT | awk '{print $2}'`
has_dos=`grep -En '(^|[[:space:]])out_dos($|[[:space:]])' INPUT | awk '{print $2}'`
has_hs=`grep -En '(^|[[:space:]])out_hs($|[[:space:]])' INPUT | awk '{print $2}'`
has_r=`grep -En '(^|[[:space:]])out_r($|[[:space:]])' INPUT | awk '{print $2}'`
out_descriptor=`grep out_descriptor INPUT | awk '{print $2}' | sed s/[[:space:]]//g`
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
if [ $calculation != "nscf" ]; then
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
        total_h=`sum_file OUT.autotest/data-H`
        echo "totalHmatrix $total_h" >>$1
	total_s=`sum_file OUT.autotest/data-S`
	echo "totalSmatrix $total_s" >>$1
fi

#echo $total_band
ttot=`grep $word $running_path | awk '{print $3}'`
echo "totaltimeref $ttot" >>$1
