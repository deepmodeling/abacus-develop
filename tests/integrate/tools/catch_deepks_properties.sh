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

get_input_key_value(){
	key=$1
	inputf=$2
	value=$(awk -v key=$key '{if($1==key) a=$2} END {print a}' $inputf)
	echo $value
}

has_force=$(get_input_key_value "cal_force" "INPUT")
has_stress=$(get_input_key_value "cal_stress" "INPUT")
deepks_out_labels=$(get_input_key_value "deepks_out_labels" "INPUT")
deepks_scf=$(get_input_key_value "deepks_scf" "INPUT")
deepks_bandgap=$(get_input_key_value "deepks_bandgap" "INPUT")
deepks_v_delta=$(get_input_key_value "deepks_v_delta" "INPUT")

#--------------------------------------------
# deepks
#--------------------------------------------
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
if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 2 ]; then
	deepks_atom=`python3 get_sum_numpy.py OUT.autotest/deepks_atom.npy `
	echo "deepks_atom $deepks_atom" >> $1
	deepks_box=`python3 get_sum_numpy.py OUT.autotest/deepks_box.npy `
	echo "deepks_box $deepks_box" >> $1
	deepks_energy=`python3 get_sum_numpy.py OUT.autotest/deepks_energy.npy `
	echo "deepks_energy $deepks_energy" >> $1
	if ! test -z "$has_force" && [ $has_force == 1 ]; then
	    deepks_force=`python3 get_sum_numpy.py OUT.autotest/deepks_force.npy `
		echo "deepks_force $deepks_force" >> $1
	fi
	if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
	    deepks_stress=`python3 get_sum_numpy.py OUT.autotest/deepks_stress.npy `
		echo "deepks_stress $deepks_stress" >> $1
	fi
	if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ]; then
	    deepks_orbital=`python3 get_sum_numpy.py OUT.autotest/deepks_orbital.npy `
		echo "deepks_orbital $deepks_orbital" >> $1
	fi
	if ! test -z "$deepks_v_delta" && [[ $deepks_v_delta -gt 0 ]]; then
	    deepks_hamiltonian=`python3 get_sum_numpy.py OUT.autotest/deepks_hamiltonian.npy `
		echo "deepks_hamiltonian $deepks_hamiltonian" >> $1
		deepks_overlap=`python3 get_sum_numpy.py OUT.autotest/deepks_overlap.npy `
		echo "deepks_overlap $deepks_overlap" >> $1
	fi
fi

#--------------------------------------------
# band gap information
#--------------------------------------------
if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ] && [ $deepks_out_labels == 1 ]; then
	odelta=`python3 get_odelta.py`
	echo "odelta $odelta" >>$1
	oprec=`python3 get_oprec.py`
	echo "oprec $oprec" >> $1
fi


#--------------------------------------------
# check vdelta in deepks
#--------------------------------------------
if ! test -z "$deepks_v_delta" && [ $deepks_v_delta == 1 ] && [ $deepks_out_labels == 1 ]; then
	totalh=`python3 get_sum_numpy.py OUT.autotest/deepks_htot.npy `
	echo "totalh $totalh" >>$1
	totalvdelta=`python3 get_v_delta.py`
	echo "totalvdelta $totalvdelta" >>$1
	totalvdp=`python3 get_sum_numpy.py OUT.autotest/deepks_vdpre.npy `
	echo "totalvdp $totalvdp" >> $1
fi

#--------------------------------------------
# check vdelta in deepks
#--------------------------------------------
if ! test -z "$deepks_v_delta" && [ $deepks_v_delta == 2 ] && [ $deepks_out_labels == 1 ]; then
	totalh=`python3 get_sum_numpy.py OUT.autotest/deepks_htot.npy `
	echo "totalh $totalh" >>$1
	totalvdelta=`python3 get_v_delta.py`
	echo "totalvdelta $totalvdelta" >>$1
	total_phialpha=`python3 get_sum_numpy.py OUT.autotest/deepks_phialpha.npy `
	echo "total_phialpha $total_phialpha" >> $1
	total_gevdm=`python3 get_sum_numpy.py OUT.autotest/deepks_gevdm.npy `
	echo "total_gevdm $total_gevdm" >> $1
fi