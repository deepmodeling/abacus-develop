#!/bin/bash

# Common library for the integrate-test property collectors.
#
# This file is meant to be *sourced*, not executed. It provides:
#   - absolute paths to the helper tools (CompareFile.py, cube_tool.py, ...),
#     resolved from this file's own location so callers no longer depend on
#     the hard-coded "../../integrate/validation_tools" relative path;
#   - the shared helper functions (sum_file, get_input_key_value,
#     sanitize_result_key, record_compare_result);
#   - props_init(): one-time parsing of all INPUT switch keys into global
#     shell variables, plus clearing of the result file;
#   - props_finalize(): the always-emitted total-time entry.
#
# Extension contract: the entry script (catch_properties.sh) sources this
# file and then each props_<category>.sh module. A module defines
# run_<category>_props() that reads the global switch variables set by
# props_init() and appends "key value" lines to the result file. Keys used
# by a single module only (e.g. deepks_*) may be read by that module itself.

# Absolute path of the directory containing this file
# (i.e. integrate/validation_tools).
PROPS_TOOLS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

COMPARE_SCRIPT="$PROPS_TOOLS_DIR/CompareFile.py"
CUBE_TOOL="$PROPS_TOOLS_DIR/cube_tool.py"
COLLECT_NPY_MEANS="$PROPS_TOOLS_DIR/collect_npy_means.py"


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


sanitize_result_key(){
	echo "$1" | sed 's/[^A-Za-z0-9_]/_/g'
}


record_compare_result(){
	result_file=$1
	result_key=$2
	ref_file=$3
	cal_file=$4
	accuracy=${5:-8}
	use_abs=${6:-0}

	if [ ! -f "$ref_file" ] || [ ! -f "$cal_file" ]; then
		echo "$result_key 1" >> "$result_file"
		return
	fi

	python3 $COMPARE_SCRIPT "$ref_file" "$cal_file" "$accuracy" -abs "$use_abs"
	echo "$result_key $?" >> "$result_file"
}


# props_init RESULT_FILE
#
# Reads the shared INPUT switch keys once into global variables, locates the
# running log, derives shared values (natom, is_lr, ...), and truncates the
# result file. Must be called from the test-case directory (where INPUT and
# OUT.autotest live) before any run_<category>_props function.
props_init(){
	props_result_file=$1

	# the command will ignore lines starting with #
	calculation=`grep calculation INPUT | grep -v '^#' | awk '{print $2}' | sed s/[[:space:]]//g`

	running_path=$(ls OUT.autotest/running_${calculation}*.log 2>/dev/null | head -1)
	if [ -z "$running_path" ]; then
	    echo "Error: No running log file found for calculation=$calculation in OUT.autotest/"
	    exit 1
	fi

	natom=`grep -En '(^|[[:space:]])TOTAL ATOM NUMBER($|[[:space:]])' $running_path | tail -1 | awk '{print $6}'`
	has_force=$(get_input_key_value "cal_force" "INPUT")
	has_stress=$(get_input_key_value "cal_stress" "INPUT")
	has_band=$(get_input_key_value "out_band" "INPUT")
	has_dos=$(get_input_key_value "out_dos" "INPUT")
	has_cond=$(get_input_key_value "cal_cond" "INPUT")
	out_hsk=$(get_input_key_value "out_hsk" "INPUT")
	out_hsr=$(get_input_key_value "out_hsr" "INPUT")
	has_hs=$(get_input_key_value "out_mat_hs" "INPUT")
	has_hs2=$(get_input_key_value "out_mat_hs2" "INPUT")
	out_hr_npz=$(get_input_key_value "out_hr_npz" "INPUT")
	out_hsr_npz=$(get_input_key_value "out_hsr_npz" "INPUT")
	out_dm_npz=$(get_input_key_value "out_dm_npz" "INPUT")
	if ! test -z "$out_hsk"; then
	    has_hs=$out_hsk
	fi
	if ! test -z "$out_hsr"; then
	    has_hs2=$out_hsr
	fi
	has_xc=$(get_input_key_value "out_mat_xc" "INPUT")
	has_xc2=$(get_input_key_value "out_mat_xc2" "INPUT")
	has_eband_separate=$(get_input_key_value "out_eband_terms" "INPUT")
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
	out_pband=$(get_input_key_value "out_proj_band" "INPUT")
	toW90=$(get_input_key_value "towannier90" "INPUT")
	has_mat_r=$(get_input_key_value "out_mat_r" "INPUT")
	has_mat_t=$(get_input_key_value "out_mat_t" "INPUT")
	has_mat_syns=$(get_input_key_value "cal_syns" "INPUT")
	has_mat_dh=$(get_input_key_value "out_mat_dh" "INPUT")
	has_scan=$(get_input_key_value "dft_functional" "INPUT")
	out_chg=$(get_input_key_value "out_chg" "INPUT")
	has_ldos=$(get_input_key_value "out_ldos" "INPUT")
	esolver_type=$(get_input_key_value "esolver_type" "INPUT")
	rdmft=$(get_input_key_value "rdmft" "INPUT")
	word_total_time="atomic_world"
	symmetry=$(get_input_key_value "symmetry" "INPUT")
	out_current=$(get_input_key_value "out_current" "INPUT")
	out_efield=$(get_input_key_value "out_efield" "INPUT")
	out_vecpot=$(get_input_key_value "out_vecpot" "INPUT")
	nspin=$(get_input_key_value "nspin" "INPUT")
	out_alllog=$(get_input_key_value "out_alllog" "INPUT")
	test -e $props_result_file && rm $props_result_file

	# if NOT non-self-consistent calculations or linear response
	is_lr=0
	if [ ! -z $esolver_type ] && ([ $esolver_type == "lr" ] || [ $esolver_type == "ks-lr" ]); then
		is_lr=1
	fi
}


# props_finalize RESULT_FILE
#
# Appends the always-present total-time entry. Kept as the single trailing
# step so new modules can be inserted before it without touching the ending.
props_finalize(){
	props_result_file=$1
	ttot=`grep $word_total_time $running_path | awk '{print $3}'`
	echo "totaltimeref $ttot" >> $props_result_file
}
