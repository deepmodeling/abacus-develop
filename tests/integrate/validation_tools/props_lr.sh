#!/bin/bash

# LR-TDDFT collectors, sourced after props_common.sh and props_init().
# Tracking cases emit one invariant check and bypass single-point collectors.
# Return success only when a tracking reference was handled.
run_lr_tracking_props(){
    if [ "$is_lr" != 1 ] || [ ! -f state_tracking.ref ]; then
        return 1
    fi
    if python3 "$PROPS_SCRIPT_DIR/excited_state_check.py" "$running_path" state_tracking.ref; then
        echo "CompareStateTracking_pass 0" >> "$props_result_file"
    else
        echo "CompareStateTracking_pass 1" >> "$props_result_file"
    fi
    return 0
}

# Reuse the generic force property; LR cases set force_threshold in their
# per-case threshold file, independently of excitation-energy accuracy.
run_lr_force_props(){
# LR prints a separate force table for each excited state and spin channel.
if [ $is_lr == 1 ] && ! test -z "$has_force" && [ $has_force == 1 ]; then
	awk '/Forces \(-gradients\) of each excited state/{flag=1; next} flag{for(i=1;i<=NF;i++) if($i=="force"){print $(i+1),$(i+2),$(i+3)}}' "$running_path" > lr_force.txt
	# Accumulate all components before rounding the LR force total.
	total_lr_force=$(awk '{for (i=1; i<=NF; ++i) sum += sqrt($i*$i)} END {printf "%.10f\n", sum}' lr_force.txt)
	rm lr_force.txt
	echo "totalforceref $total_lr_force" >>$props_result_file
	# Keep the historical total, and independently check signed, labelled components.
	local force_tolerance=0.000001
	if [ -f threshold ]; then
		local configured_force_tolerance
		configured_force_tolerance=$(get_input_key_value force_threshold threshold)
		if [ -n "$configured_force_tolerance" ]; then
			force_tolerance=$configured_force_tolerance
		fi
	fi
	if python3 "$PROPS_SCRIPT_DIR/lr_force_check.py" "$running_path" lr_force.ref "$force_tolerance"; then
		echo "CompareLRForce_pass 0" >> "$props_result_file"
	else
		echo "CompareLRForce_pass 1" >> "$props_result_file"
	fi
fi

}
