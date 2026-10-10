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

run_basic_props(){

#----------------------------
# total energy information
#----------------------------
if [ $calculation != "get_wf" ]\
&& [ $calculation != "get_pchg" ] && [ $calculation != "get_s" ]\
&& [ $is_lr == 0 ]; then
	etot=$(grep "ETOT_" "$running_path" | tail -1 | awk '{print $2}')
	etotperatom=`awk 'BEGIN {x='$etot';y='$natom';printf "%.10f\n",x/y}'`
	echo "etotref $etot" >>$props_result_file
	echo "etotperatomref $etotperatom" >>$props_result_file
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
if ! test -z "$has_force" && [ $has_force == 1 ]; then
	nn3=`echo "$natom + 3" |bc`
    # check the last step result
    grep -A$nn3 "TOTAL-FORCE" $running_path |awk 'NF==4{print $2,$3,$4}' | tail -$natom > force.txt
	total_force=`sum_file force.txt`
    rm force.txt
	echo "totalforceref $total_force" >>$props_result_file
fi

#-------------------------------
# stress information
#-------------------------------
if ! test -z "$has_stress" && [  $has_stress == 1 ]; then
    grep -A6 "TOTAL-STRESS" $running_path| awk 'NF==3' | tail -3> stress.txt
	total_stress=`sum_file stress.txt`
	rm stress.txt
	echo "totalstressref $total_stress" >>$props_result_file
fi

#-------------------------------
# DOS information
#-------------------------------
if ! test -z "$has_dos"  && [  $has_dos == 1 ]; then
	total_dos=`cat OUT.autotest/dos*.txt | awk 'END {print}' | awk '{print $3}'`
	echo "totaldosref $total_dos" >> $props_result_file
fi

#-------------------------------
# Onsager coefficiency
#-------------------------------
if ! test -z "$has_cond"  && [  $has_cond == 1 ]; then
	onref=refOnsager.txt
	oncal=OUT.autotest/Onsager.txt
	python3 $COMPARE_SCRIPT $onref $oncal 3 -com_type 0
    echo "CompareH_Failed $?" >>$props_result_file
	rm -f je-je.txt Chebycoef
fi

}

run_basic_props_post_ml(){

#--------------------------------------------
# implicit solvation model
#--------------------------------------------
if ! test -z "$imp_sol" && [ $imp_sol == 1 ]; then
	esol_el=`grep E_sol_el $running_path | awk '{print $3}'`
	esol_cav=`grep E_sol_cav $running_path | awk '{print $3}'`
	echo "esolelref $esol_el" >>$props_result_file
	echo "esolcavref $esol_cav" >>$props_result_file
fi

}

run_basic_props_post_deepks(){

#--------------------------------------------
# check symmetry
#--------------------------------------------
if ! test -z "$symmetry" && [ $symmetry == 1 ]; then
	# exclude the nspin=4 MAGNETIC POINT/SPACE GROUP lines so they do not interfere
	# with the crystallographic point-group / space-group detection below
	pointgroup=`grep 'POINT GROUP =' $running_path | grep -v 'MAGNETIC' | grep -v 'BvK' | awk '{print $4}'`
	spacegroup=`grep 'SPACE GROUP =' $running_path | grep -v 'MAGNETIC' | grep -v 'BvK' | awk '{print $7}'`
	nksibz=`grep 'Number of irreducible k-points' $running_path | awk '{print $6}'`
	echo "pointgroupref $pointgroup" >>$props_result_file
	echo "spacegroupref $spacegroup" >>$props_result_file
	echo "nksibzref $nksibz" >>$props_result_file
	# (nspin=4) magnetic (Shubnikov) group analysis: capture the space-group-consistent
	# magnetic point group. Only printed when the group is actually reduced (magnetic);
	# non-magnetic nspin=4 does not print it, so the capture is skipped when empty.
	if ! test -z "$nspin" && [ $nspin == 4 ]; then
		magpointgroup=`grep 'MAGNETIC POINT GROUP IN SPACE GROUP' $running_path | awk '{print $NF}'`
		if ! test -z "$magpointgroup"; then
			echo "magpointgroupref $magpointgroup" >>$props_result_file
		fi
	fi
fi

}

run_basic_props_post_rdmft(){

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
    log_files=$(ls OUT.autotest/running*.log 2>/dev/null)

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
        echo "log_filename_validation 1" >>$props_result_file
    else
        echo "Error: Some log filenames do not contain 'running_${calculation}_'"
        echo "log_filename_validation 0" >>$props_result_file
        exit 1
    fi
fi

}
