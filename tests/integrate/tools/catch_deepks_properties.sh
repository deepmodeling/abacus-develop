#!/bin/bash
COMPARE_SCRIPT="../../integrate/tools/CompareFile.py"

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

# General function to process npy files
process_npy() {
    local mode=$1
    local op=$2
    local file_prefix=$3
    local base_prefix=$4
    local output_name=$5
    
    local total="0"
    local file_pattern=""
    local base_pattern=""
    
    # Determine file pattern based on mode
    if [ "$mode" = "multi" ]; then
        # multi mode: multiple files with _e* suffix from different electronic steps
        file_pattern="OUT.autotest/DeePKS_Labels_Elec/${file_prefix}_e*.npy"
        base_pattern="OUT.autotest/DeePKS_Labels_Elec/${base_prefix}_e*.npy"
    elif [ "$mode" = "single" ];then
        # single mode: single file
        file_pattern="OUT.autotest/deepks_${file_prefix}.npy"
        base_pattern="OUT.autotest/deepks_${base_prefix}.npy"
    fi
    
    # Process files
    for file in $file_pattern; do
        if [ ! -f "$file" ]; then
            continue
        fi
        
        # Get step number for multi mode
        local step=""
        if [ "$mode" = "multi" ]; then
            step=$(basename "$file" | grep -oP 'e\d+')
        fi
        
        # Get corresponding base file
        local base_file=""
        if [ "$mode" = "multi" ]; then
            base_file="OUT.autotest/DeePKS_Labels_Elec/${base_prefix}_${step}.npy"
        else
            base_file=$base_pattern
        fi
        
        # Calculate value based on operation
        if [ "$op" = "abs" ]; then
            val=$(python3 ../tools/get_sum_abs.py "$file")
        elif [ "$op" = "delta" ] && [ -f "$base_file" ]; then
            val=$(python3 ../tools/get_sum_delta.py "$file" "$base_file")
        elif [ "$op" = "numpy" ]; then
            val=$(python3 ../tools/get_sum_numpy.py "$file")
        else
            val=0
        fi
        
		if [ "$mode" = "multi" ]; then
			#echo "total: $total, val: $val"
         	total=$(echo "$total + $val" | bc)
		else
			total=$val
		fi
    done
    
    echo "$output_name $total" >>$6
}

# Process a group of label outputs
process_many_npys() {
    local mode=$1
    local suffix=$2
    local output_file=$3
    
    # energy 
    process_npy "$mode" "abs" "etot" "" "deepks_e_label$suffix" "$output_file"
    process_npy "$mode" "delta" "etot" "ebase" "deepks_edelta$suffix" "$output_file"
    
    # For deepks_bandgap = 1
    if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ]; then
        process_npy "$mode" "abs" "otot" "" "deepks_o_label$suffix" "$output_file"
        process_npy "$mode" "delta" "otot" "obase" "deepks_odelta$suffix" "$output_file"
        process_npy "$mode" "numpy" "orbpre" "" "deepks_oprec$suffix" "$output_file"
    fi
    
    # For deepks_v_delta > 0
    if ! test -z "$deepks_v_delta" && [ $deepks_v_delta -gt 0 ]; then
        process_npy "$mode" "abs" "htot" "" "deepks_h_label$suffix" "$output_file"
        process_npy "$mode" "delta" "htot" "hbase" "deepks_vdelta$suffix" "$output_file"
        
        if [ $deepks_v_delta == 1 ]; then
            process_npy "$mode" "abs" "vdpre" "" "deepks_vdp$suffix" "$output_file"
        elif [ $deepks_v_delta == 2 ]; then
            process_npy "$mode" "abs" "phialpha" "" "deepks_phialpha$suffix" "$output_file"
            process_npy "$mode" "numpy" "gevdm" "" "deepks_gevdm$suffix" "$output_file"
        fi
    fi
}

# Main script
has_force=$(get_input_key_value "cal_force" "INPUT")
has_stress=$(get_input_key_value "cal_stress" "INPUT")
deepks_out_labels=$(get_input_key_value "deepks_out_labels" "INPUT")
deepks_scf=$(get_input_key_value "deepks_scf" "INPUT")
deepks_bandgap=$(get_input_key_value "deepks_bandgap" "INPUT")
deepks_v_delta=$(get_input_key_value "deepks_v_delta" "INPUT")
deepks_out_freq_elec=$(get_input_key_value "deepks_out_freq_elec" "INPUT")

#---------------------------------------------------------------------------
# Test for descriptor
#---------------------------------------------------------------------------
if ! test -z "$deepks_scf" && [ $deepks_scf == 1 ]; then
    # Process descriptor data
    sed '/n_des/d' OUT.autotest/deepks_desc.dat > des_tmp.txt
    total_des=$(sum_file des_tmp.txt 5)
    rm des_tmp.txt
    echo "deepks_desc $total_des" >>$1
    
    process_npy "single" "abs" "dm_eig" "" "deepks_dm_eig" "$1"
fi

#---------------------------------------------------------------------------
# Test for deepks_out_labels = 1 and deepks_out_freq_elec > 0
#---------------------------------------------------------------------------
if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 1 ]; then
    process_many_npys "single" "" "$1"

	# force and stress not considered in deepks_out_freq_elec > 0, so not in process_many_npys
    # For cal_force = 1
    if ! test -z "$has_force" && [ $has_force == 1 ]; then
        process_npy "single" "abs" "ftot" "" "deepks_f_label" "$1"
        process_npy "single" "delta" "ftot" "fbase" "deepks_fdelta" "$1"
        process_npy "single" "abs" "gradvx" "" "deepks_fpre" "$1"
    fi
    
    # For cal_stress = 1
    if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
        process_npy "single" "abs" "stot" "" "deepks_s_label" "$1"
        process_npy "single" "delta" "stot" "sbase" "deepks_sdelta" "$1"
        process_npy "single" "abs" "gvepsl" "" "deepks_spre" "$1"
    fi

    # For deepks_v_delta < 0
    if ! test -z "$deepks_v_delta" && [ $deepks_v_delta -lt 0 ]; then
        python3 $COMPARE_SCRIPT "deepks_hrtot.csr.ref" "OUT.autotest/deepks_hrtot.csr" 8
        echo "deepks_hr_label_pass $?" >>$1
        python3 $COMPARE_SCRIPT "deepks_hrdelta.csr.ref" "OUT.autotest/deepks_hrdelta.csr" 8
        echo "deepks_vdelta_r_pass $?" >>$1
        # For deepks_v_delta = -1
        if [ $deepks_v_delta -eq -1 ]; then
            process_npy "single" "abs" "vdrpre" "" "deepks_vdrp" "$1"
        fi
        # For deepks_v_delta = -2
        if [ $deepks_v_delta -eq -2 ]; then
            process_npy "single" "abs" "phialpha_r" "" "deepks_phialpha_r" "$1"
            process_npy "single" "numpy" "gevdm" "" "deepks_gevdm" "$1"
        fi
    fi
    
    # Process deepks_out_freq_elec > 0
    if [ ! -z "$deepks_out_freq_elec" ] && [ $deepks_out_freq_elec -gt 0 ]; then
        process_many_npys "multi" "_elec" "$1"
    fi
fi

#---------------------------------------------------------------------------
# Test for deepks_out_labels = 2
#---------------------------------------------------------------------------
if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 2 ]; then
    process_npy "single" "numpy" "atom" "" "deepks_atom" "$1"
    process_npy "single" "numpy" "box" "" "deepks_box" "$1"
    process_npy "single" "numpy" "energy" "" "deepks_energy" "$1"
    
    if ! test -z "$has_force" && [ $has_force == 1 ]; then
        process_npy "single" "numpy" "force" "" "deepks_force" "$1"
    fi
    
    if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
        process_npy "single" "numpy" "stress" "" "deepks_stress" "$1"
    fi
    
    if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ]; then
        process_npy "single" "numpy" "orbital" "" "deepks_orbital" "$1"
    fi
    
    if ! test -z "$deepks_v_delta" && [[ $deepks_v_delta -gt 0 ]]; then
        process_npy "single" "numpy" "hamiltonian" "" "deepks_hamiltonian" "$1"
        process_npy "single" "numpy" "overlap" "" "deepks_overlap" "$1"
    fi
fi