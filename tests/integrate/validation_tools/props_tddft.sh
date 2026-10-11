#!/bin/bash

# Real-time TDDFT property collectors.
#
# Sourced by catch_properties.sh after props_common.sh; relies on the global
# switch variables set by props_init() and appends result lines to
# $props_result_file.

run_tddft_props(){

#--------------------------------------------
# check currents in rt-TDDFT
#--------------------------------------------
if ! test -z "$out_current" && [ $out_current ]; then
	current1ref=current_tot.txt.ref
	current1cal=OUT.autotest/current_tot.txt
	python3 $COMPARE_SCRIPT $current1ref $current1cal 10
	echo "CompareCurrent_pass $?" >>$props_result_file
fi

#--------------------------------------------
# Check electric fields in rt-TDDFT
#--------------------------------------------
if ! test -z "$out_efield" && [ "$out_efield" == 1 ]; then
	efield_refs=(efield_*.txt.ref)
	if [ ! -e "${efield_refs[0]}" ]; then
		echo "CompareEfieldReference_pass 1" >>$props_result_file
	else
		for efield_ref in "${efield_refs[@]}"; do
			efield_name=${efield_ref%.ref}
			efield_key=$(sanitize_result_key "$efield_name")
			record_compare_result "$props_result_file" "Compare${efield_key}_pass" "$efield_ref" "OUT.autotest/$efield_name" 8
		done
	fi
fi

#--------------------------------------------
# Check vector potential in rt-TDDFT
#--------------------------------------------
if ! test -z "$out_vecpot" && [ "$out_vecpot" == 1 ]; then
	record_compare_result "$props_result_file" "CompareVectorPot_pass" "vector_pot.txt.ref" "OUT.autotest/vector_pot.txt" 8
fi

}
