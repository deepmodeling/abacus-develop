#!/bin/bash

# Real-time TDDFT property collectors.
#
# Sourced by catch_properties.sh after props_common.sh; relies on the global
# switch variables set by props_init() and appends result lines to
# $props_result_file.

run_tddft_props() {

  #--------------------------------------------
  # Check spin observables in rt-TDDFT
  #--------------------------------------------
  for spin_file in mag_tot.txt mag_atom.txt torque_tot.txt torque_atom.txt spincur_tot.txt; do
    if [ -f "$spin_file.ref" ]; then
      spin_key=$(sanitize_result_key "$spin_file")
      record_compare_result "$props_result_file" "Compare${spin_key}_pass" "$spin_file.ref" "OUT.autotest/$spin_file" 8
    fi
  done

  #--------------------------------------------
  # Check currents in rt-TDDFT
  #--------------------------------------------
  if ! test -z "$out_current" && [ $out_current ]; then
    current1ref=current_tot.txt.ref
    current1cal=OUT.autotest/current_tot.txt
    if [ "$out_current" == 2 ]; then
      current1cal=OUT.autotest/current_tot_comm.txt
    fi
    local current_tol=1e-10
    if [ -f threshold ]; then
      local configured_tol
      configured_tol=$(get_input_key_value current_threshold threshold)
      if [ -n "$configured_tol" ]; then
        current_tol=$configured_tol
      fi
    fi
    local current_status=1
    if [ -f "$current1ref" ] && [ -f "$current1cal" ]; then
      local reference_rows
      local current_rows
      reference_rows=$(awk 'NF { count++ } END { print count+0 }' "$current1ref")
      current_rows=$(awk 'NF { count++ } END { print count+0 }' "$current1cal")
      if [ "$reference_rows" -gt 0 ] && [ "$current_rows" -ge "$reference_rows" ]; then
        python3 "$COMPARE_SCRIPT" "$current1ref" "$current1cal" -tol "$current_tol"
        current_status=$?
      fi
    fi
    echo "CompareCurrent_pass $current_status" >>$props_result_file
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
