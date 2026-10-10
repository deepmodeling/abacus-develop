#!/bin/bash

# Machine-learning / advanced-method property collectors.
#
# Sourced by catch_properties.sh after props_common.sh; relies on the global
# switch variables set by props_init() and appends result lines to
# $props_result_file.
#
# Hooks (called by the entry point at the matching positions):
#   run_ml_descriptor_props() - MLKEDF .npy descriptor means
#   run_ml_rpa_props()        - RPA total energy and librpa output compare
#   run_ml_lr_props()         - linear-response excitation energies
#   run_ml_rdmft_props()      - RDMFT energy-term extraction

run_ml_descriptor_props(){

#--------------------------------------------
# ML gene data descriptors (.npy)
#--------------------------------------------
descriptor_dir="OUT.autotest/MLKEDF_Descriptors"
if [ -d "$descriptor_dir" ]; then
	python3 $COLLECT_NPY_MEANS "$descriptor_dir" >> "$props_result_file"
fi

}

run_ml_rpa_props(){

#--------------------------------------------
# random phase approximation
#--------------------------------------------
if ! test -z "$run_rpa" && [ $run_rpa == 1 ]; then
	Etot_without_rpa=`grep Etot_without_rpa log.txt | awk 'BEGIN{FS=":"} {print $2}' `
	echo "Etot_without_rpa $Etot_without_rpa" >> $props_result_file
	rpa_outdir=$(get_input_key_value "rpa_outdir" "INPUT")
	if [ -z "$rpa_outdir" ]; then
		rpa_outdir="./OUT.librpa"
	fi
	rpa_outdir=${rpa_outdir%/}
	shopt -s nullglob
	rpa_ref_files=(refcoulomb_*.txt refCs_*.txt refshrink_sinvS_*.txt)
	if [ ${#rpa_ref_files[@]} -gt 0 ]; then
		IFS=$'\n' rpa_ref_files=($(printf '%s\n' "${rpa_ref_files[@]}" | LC_ALL=C sort))
		unset IFS
		for onref in "${rpa_ref_files[@]}"; do
			oncal_name=${onref#ref}
			oncal="$rpa_outdir/$oncal_name"
			compare_key="CompareRPA_$(sanitize_result_key "$oncal_name")_pass"
			record_compare_result "$props_result_file" "$compare_key" "$onref" "$oncal" 8 1
		done
	fi
	shopt -u nullglob
fi

}

run_ml_lr_props(){

#--------------------------------------------
# Linear response function
#--------------------------------------------
if [ $is_lr == 1 ]; then
	shopt -s nullglob
	lr_files=(OUT.autotest/trans_analysis_*_tda.dat)
	if [ ${#lr_files[@]} -gt 0 ]; then
		cat "${lr_files[@]}" | awk '/Excitation Energy/{p=1; next} p && /^[[:space:]]*[0-9]+[[:space:]]/{printf "excitationenergyref%d %.6f\n", ++n, $2} /Occupied orbital/{p=0}' >>$props_result_file
	fi
	shopt -u nullglob
fi

}

run_ml_rdmft_props(){

#--------------------------------------------
# Check RDMFT method
#--------------------------------------------
if ! test -z "$rdmft" && [[ $rdmft == 1 ]]; then
	echo "" >>$props_result_file
	echo "The following energy units are in Rydberg:" >>$props_result_file

	E_TV_RDMFT=$(grep "E_TV_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_TV_RDMFT_ref $E_TV_RDMFT" >>$props_result_file

	E_hartree_RDMFT=$(grep "E_hartree_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_hartree_RDMFT_ref $E_hartree_RDMFT" >>$props_result_file

	Exc_cwp22_RDMFT=$(grep "Exc_cwp22_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "Exc_cwp22_RDMFT_ref $Exc_cwp22_RDMFT" >>$props_result_file

	E_Ewald=$(grep "E_Ewald" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_Ewald_ref $E_Ewald" >>$props_result_file

	E_entropy=$(grep "E_entropy(-TS)" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_entropy_ref $E_entropy" >>$props_result_file

	E_descf=$(grep "E_descf" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_descf_ref $E_descf" >>$props_result_file

	Etotal_RDMFT=$(grep "Etotal_RDMFT" "$running_path" | tail -1 | awk '{print $2}')
	echo "Etotal_RDMFT_ref $Etotal_RDMFT" >>$props_result_file

	Exc_ksdft=$(grep "Exc_ksdft" "$running_path" | tail -1 | awk '{print $2}')
	echo "Exc_ksdft_ref $Exc_ksdft" >>$props_result_file

	E_exx_ksdft=$(grep "E_exx_ksdft" "$running_path" | tail -1 | awk '{print $2}')
	echo "E_exx_ksdft_ref $E_exx_ksdft" >>$props_result_file

	echo "" >>$props_result_file
fi

}
