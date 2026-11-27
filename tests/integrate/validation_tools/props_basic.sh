#!/bin/bash

# Basic (energy/force/stress/DOS/Onsager/solvation/symmetry/alllog)
# property collectors.
#
# Sourced by catch_properties.sh after props_common.sh; relies on the global
# variables set by props_init() ($calculation, $running_path, $natom,
# $nspin, $is_lr, $has_*, $symmetry, $imp_sol, ...).
#
# The original monolithic script interleaved these checks with blocks that
# now live in other modules (ML descriptors, DeePKS, RDMFT). To preserve the
# exact emission order of result lines, this module exposes four hooks that
# the entry point calls at the matching positions:
#   run_basic_props()             - before the matrix/cube collectors
#   run_basic_props_post_ml()     - after the ML descriptor block
#   run_basic_props_post_deepks() - after the DeePKS block
#   run_basic_props_post_rdmft()  - after the RDMFT block
#
# All result lines are appended to $props_result_file.

run_basic_props() {

  #----------------------------
  # total energy information
  #----------------------------
  if [ $calculation != "get_wf" ] && [ $calculation != "get_pchg" ] && [ $calculation != "get_s" ] && [ $is_lr == 0 ]; then
    etot=$(grep "ETOT_" "$running_path" | tail -1 | awk '{print $2}')
    etotperatom=$(awk 'BEGIN {x='$etot';y='$natom';printf "%.10f\n",x/y}')
    echo "etotref $etot" >> $props_result_file
    echo "etotperatomref $etotperatom" >> $props_result_file
  fi

  # Opt-in collinear magnetic-state check. Comparing both moments distinguishes
  # compensated AFM from NM; the reference tolerance is 1e-3 mu_B per cell.
  if [ "$nspin" = "2" ] && [ -f magnetism.ref ]; then
    awk '
        /Total magnetism \(Bohr mag\/cell\)/ {total = $NF; have_total = 1}
        /Absolute magnetism \(Bohr mag\/cell\)/ {absolute = $NF; have_absolute = 1}
        END {
            if (!have_total || !have_absolute) exit 1
            print total, absolute
        }
    ' "$running_path" > magnetism.out
    record_compare_result "$props_result_file" "CompareMagnetism_pass" "magnetism.ref" "magnetism.out" 3
  fi

  #----------------------------
  # force information
  #----------------------------
  if [ -n "$has_force" ] && [ "$has_force" = 1 ] && [ "$is_lr" = 0 ]; then
    if ! [[ "$natom" =~ ^[0-9]+$ ]] || [ "$natom" -eq 0 ]; then
      echo "[ERROR] Invalid atom count for force collection: $natom" >&2
      return 1
    fi
    # The final MD structure contains the evaluated force on every atom.
    if [ "$calculation" = "md" ]; then
      force_source=$(awk '/Write file: .*STRU_FINAL$/ {sub(/^.*Write file: /, ""); path=$0} END {print path}' "$running_path")
      if [ -z "$force_source" ] || [ ! -f "$force_source" ]; then
        force_source=OUT.autotest/STRU_FINAL
        if [ ! -f "$force_source" ] && [ -f OUT.autotest/STRUC/STRU_FINAL ]; then
          force_source=OUT.autotest/STRUC/STRU_FINAL
        fi
      fi
      force_format=stru
    else
      force_source=$running_path
      force_format=log
    fi
    if ! awk -v natom="$natom" -v format="$force_format" '
      function number(value) {
        return value ~ /^[-+]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][-+]?[0-9]+)?$/
      }
      function reset_frame() {
        in_atoms = 0
        atom_count = 0
        invalid = 0
        delete force_x
        delete force_y
        delete force_z
      }
      format == "stru" && /^ATOMIC_POSITIONS/ { reset_frame(); in_atoms = 1; next }
      format == "log" && /TOTAL-FORCE/ { reset_frame(); in_atoms = 1; next }
      in_atoms && ((format == "stru" && NF >= 3 && number($1) && number($2) && number($3)) || (format == "log" && NF == 4 && $1 != "Atoms")) {
        x = $2
        y = $3
        z = $4
        if (format == "stru") {
          found_force = 0
          for (field = 4; field <= NF; ++field) {
            if ($field == "f") {
              found_force = 1
              x = $(field + 1)
              y = $(field + 2)
              z = $(field + 3)
              break
            }
          }
          if (!found_force) invalid = 1
        }
        if (!number(x) || !number(y) || !number(z)) invalid = 1
        atom_count++
        force_x[atom_count] = x
        force_y[atom_count] = y
        force_z[atom_count] = z
        if (atom_count == natom) in_atoms = 0
      }
      END {
        if (atom_count != natom || invalid) exit 1
        for (iat = 1; iat <= natom; ++iat) {
          print force_x[iat], force_y[iat], force_z[iat]
        }
      }
    ' "$force_source" > force.txt; then
      echo "[ERROR] Failed to extract $natom valid force rows from $force_source" >&2
      rm -f force.txt
      return 1
    fi
    total_force=$(sum_file force.txt)
    rm force.txt
    echo "totalforceref $total_force" >> $props_result_file
  fi

  #-------------------------------
  # stress information
  #-------------------------------
  if ! test -z "$has_stress" && [ $has_stress == 1 ]; then
    if [ "$calculation" = "md" ]; then
      final_structure="OUT.autotest/STRU_FINAL"
      if [ ! -f "$final_structure" ]; then
        final_structure="OUT.autotest/STRUC/STRU_FINAL"
      fi
      if [ ! -f "$final_structure" ]; then
        final_structure="OUT.autotest/STRU_FINAL.cif"
      fi
      if [ ! -f "$final_structure" ]; then
        final_structure="OUT.autotest/STRUC/STRU_FINAL.cif"
      fi
      if [ ! -f "$final_structure" ]; then
        echo "[ERROR] Cannot identify the final MD step without STRU_FINAL" >&2
        return 1
      fi
      final_step=$(awk '/^# MD step [0-9]+/ { print $4; exit }' "$final_structure")
      if [ -z "$final_step" ]; then
        echo "[ERROR] Missing MD step marker in $final_structure" >&2
        return 1
      fi
      # The MD virial is the electronic stress after external pressure is
      # subtracted. Reconstruct the solver's former TOTAL-STRESS value so
      # existing reference data continue to compare the same quantity.
      press1=$(get_input_key_value press1 INPUT)
      press2=$(get_input_key_value press2 INPUT)
      press3=$(get_input_key_value press3 INPUT)
      if [ "$esolver_type" = "ofdft" ] || [ "$esolver_type" = "sdft" ]; then
        # These solvers return their electronic stress without subtracting
        # the requested external pressure.
        press1=0
        press2=0
        press3=0
      fi
      if ! awk -v p1="$press1" -v p2="$press2" -v p3="$press3" -v final="$final_step" '
        /MD STEP:/ { step = $3 }
        /VIRIAL \(kbar\):/ { rows = 3; count = 0; selected = (step == final); next }
        rows > 0 {
          if (NF != 3) exit 1
          count++
          if (selected) value[count] = $1 " " $2 " " $3
          rows--
        }
        END {
          if (!selected || count != 3 || rows != 0) exit 1
          split(value[1], first)
          split(value[2], second)
          split(value[3], third)
          # Keep the precision of the former TOTAL-STRESS table. awk print
          # defaults to six significant digits, enough to change the summed
          # stress by more than the existing test tolerance.
          printf "%.10f %.10f %.10f\n", first[1] + p1, first[2], first[3]
          printf "%.10f %.10f %.10f\n", second[1], second[2] + p2, second[3]
          printf "%.10f %.10f %.10f\n", third[1], third[2], third[3] + p3
        }
      ' "$running_path" > stress.txt; then
        echo "[ERROR] Failed to extract the MD virial from $running_path" >&2
        rm -f stress.txt
        return 1
      fi
    else
      grep -A6 "TOTAL-STRESS" "$running_path" | awk 'NF==3' | tail -3 > stress.txt
    fi
    total_stress=$(sum_file stress.txt)
    rm stress.txt
    echo "totalstressref $total_stress" >> $props_result_file
  fi

  #-------------------------------
  # DOS information
  #-------------------------------
  if ! test -z "$has_dos" && [ $has_dos == 1 ]; then
    total_dos=$(cat OUT.autotest/dos*.txt | awk 'END {print}' | awk '{print $3}')
    echo "totaldosref $total_dos" >> $props_result_file
  elif [ "$has_dos" = 2 ]; then
    local dos_dir=OUT.autotest
    if [ "$out_app_flag" = 2 ]; then dos_dir+=/DOS; fi
    python3 "$PROPS_TOOLS_DIR/check_dos_npy.py" "$dos_dir" "${nspin:-1}" dos.ref >> "$props_result_file" || return 1
  fi

  #-------------------------------
  # Onsager coefficiency
  #-------------------------------
  if ! test -z "$has_cond" && [ $has_cond == 1 ]; then
    onref=refOnsager.txt
    oncal=OUT.autotest/Onsager.txt
    python3 $COMPARE_SCRIPT $onref $oncal 3 -com_type 0
    echo "CompareH_Failed $?" >> $props_result_file
    rm -f je-je.txt Chebycoef
  fi

}

run_basic_props_post_ml() {

  #--------------------------------------------
  # implicit solvation model
  #--------------------------------------------
  if ! test -z "$imp_sol" && { [ "$imp_sol" = 1 ] || [ "$imp_sol" = 2 ]; }; then
    esol_el=$(grep E_sol_el $running_path | tail -1 | awk '{print $3}')
    esol_cav=$(grep E_sol_cav $running_path | tail -1 | awk '{print $3}')
    echo "esolelref $esol_el" >> $props_result_file
    echo "esolcavref $esol_cav" >> $props_result_file
  fi

  # Solvent-aware SCCS (imp_sol 2, sccs_solvent_radius > 0). With sccs_debug 2
  # the running log must report the filled volume int(s_sa - s) dV of the last
  # evaluation, which is compared.
  solvent_aware=0
  if [ "$imp_sol" = 2 ] && awk -v radius="$sccs_solvent_radius" 'BEGIN {exit !(radius + 0 > 0)}'; then
    solvent_aware=1
  fi
  if [ "$solvent_aware" = 1 ] && [ "$sccs_debug" = 2 ]; then
    filled_line=$(grep "V_FILLED/Bohr3" "$running_path" | tail -1)
    filled_volume=$(echo "$filled_line" | awk '{for (i = 1; i < NF; ++i) if ($i == "V_FILLED/Bohr3") print $(i + 1)}')
    if [ -z "$filled_volume" ]; then
      echo "Error: solvent-aware SCCS with sccs_debug 2 but no V_FILLED/Bohr3 in $running_path"
      exit 1
    fi
    echo "sccs_filled_volumeref $filled_volume" >> $props_result_file
  fi

  # out_sol 1: the cubes must exist (the cases set imp_sol 2, or imp_sol 1 with
  # eb_k > 1, so sol_cavity.cube is written; the solvent-aware cavity adds two).
  # The cavity volume int s dV of sol_cavity.cube is compared; the other cubes
  # are only checked for presence: the sol_eps.cube integral is dominated by the
  # bulk permittivity times the volume, and eps follows from the cavity.
  if [[ "$out_sol" == "1" ]]; then
    sol_cubes="sol_eps.cube sol_cavity.cube"
    if [ "$solvent_aware" = 1 ]; then
      sol_cubes="$sol_cubes sol_cavity_local.cube sol_filled_fraction.cube"
    fi
    for cube in $sol_cubes; do
      if [ ! -f "OUT.autotest/$cube" ]; then
        echo "Error: out_sol 1 but OUT.autotest/$cube is missing"
        exit 1
      fi
    done
    # Append mode retains earlier geometries; compare the final complete frame.
    integral=$(python3 "$CUBE_TOOL" integrate-last "OUT.autotest/sol_cavity.cube") || return 1
    echo "sol_cavity_integralref $integral" >> $props_result_file
  fi

  #--------------------------------------------
  # vacuum point-countercharge correction
  #--------------------------------------------
  if [[ "$assume_isolated" == "pcc_0d" || "$assume_isolated" == "pcc_2d" ]]; then
    epcc=$(grep E_pcc "$running_path" | tail -1 | awk '{print $3}')
    echo "epccref $epcc" >> $props_result_file
  fi

}

run_basic_props_post_deepks() {

  #--------------------------------------------
  # check symmetry
  #--------------------------------------------
  if ! test -z "$symmetry" && [ $symmetry == 1 ]; then
    # exclude the nspin=4 MAGNETIC POINT/SPACE GROUP lines so they do not interfere
    # with the crystallographic point-group / space-group detection below
    pointgroup=$(grep 'POINT GROUP =' $running_path | grep -v 'MAGNETIC' | grep -v 'BvK' | awk '{print $4}')
    spacegroup=$(grep 'SPACE GROUP =' $running_path | grep -v 'MAGNETIC' | grep -v 'BvK' | awk '{print $7}')
    nksibz=$(grep 'Number of irreducible k-points' $running_path | awk '{print $6}')
    echo "pointgroupref $pointgroup" >> $props_result_file
    echo "spacegroupref $spacegroup" >> $props_result_file
    echo "nksibzref $nksibz" >> $props_result_file
    # (nspin=4) magnetic (Shubnikov) group analysis: capture the space-group-consistent
    # magnetic point group. Only printed when the group is actually reduced (magnetic);
    # non-magnetic nspin=4 does not print it, so the capture is skipped when empty.
    if ! test -z "$nspin" && [ $nspin == 4 ]; then
      magpointgroup=$(grep 'MAGNETIC POINT GROUP IN SPACE GROUP' $running_path | awk '{print $NF}')
      if ! test -z "$magpointgroup"; then
        echo "magpointgroupref $magpointgroup" >> $props_result_file
      fi
    fi
  fi

}

run_basic_props_post_rdmft() {

  #--------------------------------------------
  # Check if out_alllog is set to 1
  # and verify running*.log filenames
  #--------------------------------------------
  if ! test -z "$out_alllog" && [ $out_alllog -eq 1 ]; then
    if [ -z "$calculation" ]; then
      echo "Error: calculation parameter not found in INPUT"
      exit 1
    fi

    # Find all running*.log files in OUT.autotest directory
    log_files=$(ls OUT.autotest/running*.log 2> /dev/null)

    if [ -z "$log_files" ]; then
      echo "Error: No running*.log files found in OUT.autotest/"
      exit 1
    fi

    # Check each log file name contains the calculation parameter
    all_valid=true
    for log_file in $log_files; do
      filename=$(basename "$log_file")
      if [[ ! "$filename" =~ running_${calculation}_ ]]; then
        echo "Error: Invalid log filename $filename - should contain 'running_${calculation}_'"
        all_valid=false
      fi
    done

    if $all_valid; then
      echo "All log filenames contain 'running_${calculation}_' - validation passed"
      echo "log_filename_validation 1" >> $props_result_file
    else
      echo "Error: Some log filenames do not contain 'running_${calculation}_'"
      echo "log_filename_validation 0" >> $props_result_file
      exit 1
    fi
  fi

}
