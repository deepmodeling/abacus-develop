#!/bin/bash

# LR-TDDFT collectors, sourced after props_common.sh and props_init().
# Tracking cases emit one invariant check and bypass single-point collectors.
# Return success only when a tracking reference was handled.
run_lr_tracking_props() {
  if [ "$is_lr" != 1 ] || [ ! -f state_tracking.ref ]; then
    return 1
  fi
  run_and_record_status "$props_result_file" "CompareStateTracking_pass" \
    python3 "$PROPS_SCRIPT_DIR/excited_state_check.py" "$running_path" state_tracking.ref
  return 0
}

# Reuse the generic force property; LR cases set force_threshold in their
# per-case threshold file, independently of excitation-energy accuracy.
run_lr_force_props() {
  # LR prints a separate force table for each excited state and spin channel.
  if [ $is_lr == 1 ] && ! test -z "$has_force" && [ $has_force == 1 ]; then
    awk '
      /Forces \(-gradients\) of each excited state:|Excitation-energy force contributions \(-dOmega\/dR\) of each state:/ {
        active = 1
        next
      }
      active && (NF == 0 || ($1 == "state" && $2 == "atom")) { next }
      active && ($2 == "force" || $3 == "force") {
        for (i = 1; i <= NF; ++i) {
          if ($i == "force") print $(i+1), $(i+2), $(i+3)
        }
        next
      }
      active { active = 0 }
    ' "$running_path" > lr_force.txt
    # Accumulate all components before rounding the LR force total.
    total_lr_force=$(awk '{for (i=1; i<=NF; ++i) sum += sqrt($i*$i)} END {printf "%.10f\n", sum}' lr_force.txt)
    rm lr_force.txt
    echo "totalforceref $total_lr_force" >> $props_result_file
    # Keep the historical total, and independently check signed, labelled components.
    local force_tolerance=0.000001
    if [ -f threshold ]; then
      local configured_force_tolerance
      configured_force_tolerance=$(get_input_key_value force_threshold threshold)
      if [ -n "$configured_force_tolerance" ]; then
        force_tolerance=$configured_force_tolerance
      fi
    fi
    run_and_record_status "$props_result_file" "CompareLRForce_pass" \
      python3 "$PROPS_SCRIPT_DIR/lr_force_check.py" "$running_path" lr_force.ref "$force_tolerance"
  fi

}

# Optional exact filename check for channel-labelled electron-hole densities.
run_lr_density_props() {
  if [ "$is_lr" != 1 ] || [ ! -f lr_density.ref ]; then
    return 0
  fi
  local expected_names
  local actual_names=""
  local density_file
  local valid=1
  expected_names=$(sed '/^[[:space:]]*$/d' lr_density.ref | sort)
  for density_file in OUT.autotest/LR_e-h_density_*.cube; do
    if [ ! -s "$density_file" ]; then
      valid=0
    else
      actual_names+="${density_file##*/}"$'\n'
    fi
  done
  actual_names=$(printf '%s' "$actual_names" | sort)
  if [ "$valid" == 1 ] && [ -n "$expected_names" ] && [ "$actual_names" = "$expected_names" ]; then
    echo "CompareLRDensity_pass 0" >>"$props_result_file"
  else
    echo "LR density filenames missing, empty, or inconsistent with lr_density.ref" >&2
    echo "CompareLRDensity_pass 1" >>"$props_result_file"
  fi
}
