#!/bin/bash

# Real-time TDDFT property collectors.
#
# Sourced by catch_properties.sh after props_common.sh; relies on the global
# switch variables set by props_init() and appends result lines to
# $props_result_file.

run_tddft_props(){

# Opt-in continuity check for MD cases with one electronic step per MD step.
# md_scf.ref lists post-initial steps as: step, iterations, converged solves.
# Also require the initial ground-state SCF to converge, without fixing its
# iteration count across different numerical backends.
if [ "$calculation" = md ] && [ "$esolver_type" = tddft ] && [ -f md_scf.ref ]; then
    md_scf_status=1
    if awk '
        function emit_step() {
            if (step == 1) initial_steps++
            if (step > 0 && converged != 1) invalid = 1
            if (step > 1) {
                print step, iterations, converged
                emitted++
            }
        }
        /STEP OF MOLECULAR DYNAMICS:/ {
            emit_step()
            step = $NF
            iterations = 0
            converged = 0
        }
        /#ELEC ITER#/ { iterations = $NF }
        /#SCF IS CONVERGED#/ { converged++ }
        END {
            emit_step()
            if (emitted == 0 || initial_steps != 1 || invalid) exit 1
        }
    ' "$running_path" > md_scf.out && cmp -s md_scf.ref md_scf.out; then
        md_scf_status=0
    fi
    echo "CompareMdScf_pass $md_scf_status" >> "$props_result_file"
fi

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
