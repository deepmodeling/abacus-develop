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

#---------------------------------------------------------------------------
# Test for deepks_scf = 1
#---------------------------------------------------------------------------
if ! test -z "$deepks_scf" && [ $deepks_scf == 1 ]; then
	sed '/n_des/d' OUT.autotest/deepks_desc.dat > des_tmp.txt
	total_des=`sum_file des_tmp.txt 5`
	rm des_tmp.txt
	echo "deepks_desc $total_des" >>$1
	if ! test -z "$deepks_scf" && [ $deepks_scf == 1 ]; then
		deepks_dm_eig=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_dm_eig.npy`
	    echo "deepks_dm_eig $deepks_dm_eig" >>$1
	fi
fi

#---------------------------------------------------------------------------
# Test for deepks_out_labels = 1 (requires deepks_scf = 1)
#---------------------------------------------------------------------------
if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 1 ]; then
	if ! test -z "$deepks_scf" && [ $deepks_scf == 1 ]; then
		deepks_e_label=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_etot.npy`
	    echo "deepks_e_label $deepks_e_label" >>$1
        deepks_edelta=`python3 ../tools/get_sum_delta.py OUT.autotest/deepks_etot.npy OUT.autotest/deepks_ebase.npy`
        echo "deepks_edelta $deepks_edelta" >>$1
        # For cal_force = 1
        if ! test -z "$has_force" && [ $has_force == 1 ]; then
            deepks_f_label=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_ftot.npy`
            echo "deepks_f_label $deepks_f_label" >>$1
            deepks_fdelta=`python3 ../tools/get_sum_delta.py OUT.autotest/deepks_ftot.npy OUT.autotest/deepks_fbase.npy`
            echo "deepks_fdelta $deepks_fdelta" >>$1
            deepks_fpre=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_gradvx.npy`
            echo "deepks_fpre $deepks_fpre" >>$1
        fi
        # For cal_stress = 1
        if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
            deepks_s_label=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_stot.npy`
            echo "deepks_s_label $deepks_s_label" >>$1
            deepks_sdelta=`python3 ../tools/get_sum_delta.py OUT.autotest/deepks_stot.npy OUT.autotest/deepks_sbase.npy`
            echo "deepks_sdelta $deepks_sdelta" >>$1
            deepks_spre=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_gvepsl.npy`
            echo "deepks_spre $deepks_spre" >>$1
        fi
        # For deepks_bandgap = 1
        if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ]; then
            deepks_o_label=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_otot.npy`
            echo "deepks_o_label $deepks_o_label" >>$1
            deepks_odelta=`python3 ../tools/get_sum_delta.py OUT.autotest/deepks_otot.npy OUT.autotest/deepks_obase.npy`
            echo "deepks_odelta $deepks_odelta" >>$1
            deepks_oprec=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_orbpre.npy`
            echo "deepks_oprec $deepks_oprec" >> $1
        fi
        # For deepks_v_delta > 0
        if ! test -z "$deepks_v_delta" && [ $deepks_v_delta > 0 ]; then
            deepks_h_label=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_htot.npy`
            echo "deepks_h_label $deepks_h_label" >>$1
            deepks_vdelta=`python3 ../tools/get_sum_delta.py OUT.autotest/deepks_htot.npy OUT.autotest/deepks_hbase.npy`
            echo "deepks_vdelta $deepks_vdelta" >>$1
            # For deepks_v_delta = 1
            if [ $deepks_v_delta == 1 ]; then
                deepks_vdp=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_vdpre.npy `
                echo "deepks_vdp $deepks_vdp" >> $1
            fi
            # For deepks_v_delta = 2
            if [ $deepks_v_delta == 2 ]; then
                deepks_phialpha=`python3 ../tools/get_sum_abs.py OUT.autotest/deepks_phialpha.npy `
                echo "deepks_phialpha $deepks_phialpha" >> $1
                deepks_gevdm=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_gevdm.npy `
                echo "deepks_gevdm $deepks_gevdm" >> $1
            fi
        fi
    else
        echo "Warning: deepks_out_labels = 1 requires deepks_scf = 1"
        exit 1
    fi
fi

#---------------------------------------------------------------------------
# Test for deepks_out_labels = 2
#---------------------------------------------------------------------------
if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 2 ]; then
	deepks_atom=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_atom.npy `
	echo "deepks_atom $deepks_atom" >> $1
	deepks_box=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_box.npy `
	echo "deepks_box $deepks_box" >> $1
	deepks_energy=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_energy.npy `
	echo "deepks_energy $deepks_energy" >> $1
	if ! test -z "$has_force" && [ $has_force == 1 ]; then
	    deepks_force=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_force.npy `
		echo "deepks_force $deepks_force" >> $1
	fi
	if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
	    deepks_stress=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_stress.npy `
		echo "deepks_stress $deepks_stress" >> $1
	fi
	if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ]; then
	    deepks_orbital=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_orbital.npy `
		echo "deepks_orbital $deepks_orbital" >> $1
	fi
	if ! test -z "$deepks_v_delta" && [[ $deepks_v_delta -gt 0 ]]; then
	    deepks_hamiltonian=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_hamiltonian.npy `
		echo "deepks_hamiltonian $deepks_hamiltonian" >> $1
		deepks_overlap=`python3 ../tools/get_sum_numpy.py OUT.autotest/deepks_overlap.npy `
		echo "deepks_overlap $deepks_overlap" >> $1
	fi
fi
