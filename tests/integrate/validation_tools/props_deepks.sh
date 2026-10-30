#!/bin/bash

# DeePKS label / descriptor property collectors.
#
# Sourced by catch_properties.sh after props_common.sh; can also be sourced
# by the catch_deepks_properties.sh standalone entry. Relies on the shared
# helpers from props_common.sh ($COMPARE_SCRIPT, sum_file,
# get_input_key_value) and appends result lines to $props_result_file.
#
# The get_sum_*.py helpers below are invoked with a "../tools/" path that is
# relative to the test-case CWD (tests/09_DeePKS/<case>), where it resolves
# to tests/09_DeePKS/tools/; this is intentional and must stay CWD-relative.

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

    # Check if force or stress
    local is_force_or_stress=0
    if [ "$file_prefix" = "ftot" ] || [ "$file_prefix" = "stot" ]; then
        is_force_or_stress=1
    fi

    # Determine file pattern based on mode
    if [ "$mode" = "multi" ]; then
        # multi mode: multiple files with _e* suffix from different electronic steps, unless for force and stress
        if [ "$is_force_or_stress" -eq 1 ] ; then
            file_pattern="OUT.autotest/DeePKS_Labels_Elec/${file_prefix}.npy"
            base_pattern="OUT.autotest/DeePKS_Labels_Elec/${base_prefix}.npy"
        else
            file_pattern="OUT.autotest/DeePKS_Labels_Elec/${file_prefix}_e*.npy"
            base_pattern="OUT.autotest/DeePKS_Labels_Elec/${base_prefix}_e*.npy"
        fi
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
        # Force/stress multi-mode files (ftot.npy, stot.npy) carry no e<step>
        # suffix, so grep legitimately finds no match and returns 1. Tolerate
        # the empty match explicitly so the step extraction is also safe when
        # this module is sourced under `bash -e`.
        local step=""
        if [ "$mode" = "multi" ]; then
            step=$(basename "$file" | grep -oP 'e\d+' || true)
        fi

        # Get corresponding base file
        local base_file=""
        if [ "$mode" = "multi" ] && [ "$is_force_or_stress" -ne 1 ]; then
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

    # For deepks_bandgap > 0
    if ! test -z "$deepks_bandgap" && [ $deepks_bandgap -gt 0 ]; then
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

    if ! test -z "$has_force" && [ $has_force == 1 ]; then
        process_npy "$mode" "abs" "ftot" "" "deepks_f_label$suffix" "$output_file"
        process_npy "$mode" "delta" "ftot" "fbase" "deepks_fdelta$suffix" "$output_file"
    fi

    # For cal_stress = 1
    if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
        process_npy "$mode" "abs" "stot" "" "deepks_s_label$suffix" "$output_file"
        process_npy "$mode" "delta" "stot" "sbase" "deepks_sdelta$suffix" "$output_file"
    fi
}

run_deepks_props(){

# Main script
# cal_force/cal_stress are also parsed by props_init(); they are re-read here
# because the standalone catch_deepks_properties.sh entry skips props_init().
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
    echo "deepks_desc $total_des" >>$props_result_file

    process_npy "single" "abs" "dm_eig" "" "deepks_dm_eig" "$props_result_file"
fi

#---------------------------------------------------------------------------
# Test for deepks_out_labels = 1 and deepks_out_freq_elec > 0
#---------------------------------------------------------------------------
if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 1 ]; then
    process_many_npys "single" "" "$props_result_file"

	# gradvx and gvepsl not considered in deepks_out_freq_elec > 0, so not in process_many_npys
    # For cal_force = 1
    if ! test -z "$has_force" && [ $has_force == 1 ]; then
        process_npy "single" "abs" "gradvx" "" "deepks_fpre" "$props_result_file"
    fi

    # For cal_stress = 1
    if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
        process_npy "single" "abs" "gvepsl" "" "deepks_spre" "$props_result_file"
    fi

    # For deepks_v_delta < 0
    # CompareFile.py exits 1 when the files differ; capture the status in a
    # named variable so it is recorded even when sourced under `bash -e`.
    if ! test -z "$deepks_v_delta" && [ $deepks_v_delta -lt 0 ]; then
        hr_pass=0
        python3 $COMPARE_SCRIPT "deepks_hrtot.csr.ref" "OUT.autotest/deepks_hrtot.csr" 8 || hr_pass=$?
        echo "deepks_hr_label_pass $hr_pass" >>$props_result_file
        hrdelta_pass=0
        python3 $COMPARE_SCRIPT "deepks_hrdelta.csr.ref" "OUT.autotest/deepks_hrdelta.csr" 8 || hrdelta_pass=$?
        echo "deepks_vdelta_r_pass $hrdelta_pass" >>$props_result_file
        # For deepks_v_delta = -1
        if [ $deepks_v_delta -eq -1 ]; then
            process_npy "single" "abs" "vdrpre" "" "deepks_vdrp" "$props_result_file"
        fi
        # For deepks_v_delta = -2
        if [ $deepks_v_delta -eq -2 ]; then
            process_npy "single" "abs" "phialpha_r" "" "deepks_phialpha_r" "$props_result_file"
            process_npy "single" "numpy" "gevdm" "" "deepks_gevdm" "$props_result_file"
        fi
    fi

    # Process deepks_out_freq_elec > 0
    if [ ! -z "$deepks_out_freq_elec" ] && [ $deepks_out_freq_elec -gt 0 ]; then
        process_many_npys "multi" "_elec" "$props_result_file"
        if ! test -z "$deepks_v_delta" && [[ $deepks_v_delta -gt 0 ]]; then
            process_npy "single" "abs" "overlap" "" "deepks_overlap" "$props_result_file"
            process_npy "multi" "abs" "overlap" "" "deepks_overlap_elec" "$props_result_file"

        fi
    fi
fi

#---------------------------------------------------------------------------
# Test for deepks_out_labels = 2
#---------------------------------------------------------------------------
if ! test -z "$deepks_out_labels" && [ $deepks_out_labels == 2 ]; then
    process_npy "single" "numpy" "atom" "" "deepks_atom" "$props_result_file"
    process_npy "single" "numpy" "box" "" "deepks_box" "$props_result_file"
    process_npy "single" "numpy" "energy" "" "deepks_energy" "$props_result_file"

    if ! test -z "$has_force" && [ $has_force == 1 ]; then
        process_npy "single" "numpy" "force" "" "deepks_force" "$props_result_file"
    fi

    if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
        process_npy "single" "numpy" "stress" "" "deepks_stress" "$props_result_file"
    fi

    if ! test -z "$deepks_bandgap" && [ $deepks_bandgap == 1 ]; then
        process_npy "single" "numpy" "orbital" "" "deepks_orbital" "$props_result_file"
    fi

    if ! test -z "$deepks_v_delta" && [[ $deepks_v_delta -gt 0 ]]; then
        process_npy "single" "numpy" "hamiltonian" "" "deepks_hamiltonian" "$props_result_file"
        process_npy "single" "numpy" "overlap" "" "deepks_overlap" "$props_result_file"
    fi
fi

}
