#!/bin/bash

# Real-space cube / wave-function property collectors.
#
# Sourced by catch_properties.sh after props_common.sh; relies on the global
# switch variables set by props_init() and appends result lines to
# $props_result_file.
#
# The cube collectors originally wrapped around the out_dm matrix block
# (props_mat.sh), so this module exposes three hooks called at the matching
# positions to preserve result-line order:
#   run_cube_pot_props()  - out_pot (1/2) and out_elf cubes
#   run_cube_props()      - chg/SAN tau/LDOS and wfc real-space/PW/LCAO checks
#   run_cube_tail_props() - mulliken, pchg cubes, get_wf/get_pchg cube
#                           integration/fingerprints, nspin=4 spinor identity

run_cube_pot_props(){

#-------------------------------
# echo $out_pot1
#-------------------------------
if ! test -z "$out_pot"  && [  $out_pot == 1 ]; then
	pot1ref=pot.cube.ref
	pot1cal=OUT.autotest/pot.cube
	python3 $COMPARE_SCRIPT $pot1ref $pot1cal 3
	echo "ComparePot1_pass $?" >>$props_result_file
fi

#-------------------------------
#echo $out_pot2
#-------------------------------
if ! test -z "$out_pot"  && [  $out_pot == 2 ]; then
	pot1ref=potes.cube.ref
	pot1cal=OUT.autotest/potes.cube
	python3 $COMPARE_SCRIPT $pot1ref $pot1cal 8
	echo "ComparePot_pass $?" >>$props_result_file
fi

#-------------------------------
# Electron localized function
# echo $out_elf
#-------------------------------
if ! test -z "$out_elf"  && [  $out_elf == 1 ]; then
	elf1ref=elftot.cube.ref
	elf1cal=OUT.autotest/elftot.cube
	python3 $COMPARE_SCRIPT $elf1ref $elf1cal 3
	echo "ComparePot1_pass $?" >>$props_result_file
fi

}

run_cube_props(){

#---------------------------------------
# Charge density
#---------------------------------------
#echo $out_chg
if ! test -z "$out_chg"  && [  $out_chg -ge 1 ]; then
	record_compare_result "$props_result_file" "chg.cube_pass" "chg.cube.ref" "OUT.autotest/chg.cube" 6
fi


#---------------------------------------
# SCAN exchange-correlation information
#echo $has_scan
#---------------------------------------
if ! test -z "$has_scan"  && [  $has_scan == "scan" ] && \
       ! test -z "$out_chg" && [ $out_chg == 1 ]; then
    python3 $COMPARE_SCRIPT tau.cube.ref OUT.autotest/tau.cube 8
    echo "tau.cube_pass $?" >>$props_result_file
fi

#---------------------------------------
# local density of states
# echo $has_ldos
#---------------------------------------
if ! test -z "$has_ldos"  && [  $has_ldos == 1 ]; then
    stm_bias=$(get_input_key_value "stm_bias" "OUT.autotest/INPUT.info")
    python3 $COMPARE_SCRIPT LDOS.cube.ref OUT.autotest/LDOS_"$stm_bias"eV.cube 8
    echo "LDOS.cube_pass $?" >> $props_result_file
fi

#---------------------------------------
# wave functions in real space
# echo "$has_wfc_r" ## test out_wfc_r > 0
#---------------------------------------
if ! test -z "$has_wfc_r"  && [ $has_wfc_r == 1 ]; then
	if [[ ! -f OUT.autotest/running_scf.log ]];then
		echo "Can't find file OUT.autotest/running_scf.log"
		exit 1
	fi
	nband=$(grep NBANDS OUT.autotest/running_scf.log|awk '{print $3}')
    allgrid=$(grep "fft grid for wave functions" OUT.autotest/running_scf.log | awk -F "[=,\\\[\\\]]" '{print $3*$4*$5}')
	for((band=0;band<$nband;band++));do
		if [[ -f "OUT.autotest/wfc_realspace/wfc_realspace_0_$band" ]];then
			variance_wfc_r=`sed -n "13,$"p OUT.autotest/wfc_realspace/wfc_realspace_0_$band | \
						awk -v all=$allgrid 'BEGIN {sumall=0} {for(i=1;i<=NF;i++) {sumall+=($i-1)*($i-1)}}\
						END {printf"%.5f",(sumall/all)}'`
			echo "variance_wfc_r_0_$band $variance_wfc_r" >>$props_result_file
		else
			echo "Can't find file OUT.autotest/wfc_realspace/wfc_realspace_0_$band"
			exit 1
		fi
	done
fi

#--------------------------------------------
# wave functions in plane wave basis
# echo "$has_wfc_pw" ## test out_wfc_pw > 0
#--------------------------------------------
if ! test -z "$has_wfc_pw"  && [ $has_wfc_pw == 1 ]; then

    # according to nspin value
    if [ "$nspin" -eq 1 ]; then
        filename="wfk1_pw.txt"
    elif [ "$nspin" -eq 2 ]; then
        filename="wfk1s1_pw.txt"
    else
        # other nspin cases
        echo "Unsupported nspin value: $nspin"
        exit 1
    fi

    full_path="OUT.autotest/$filename"

	if [[ ! -f "$full_path" ]];then
		echo "Can't find file $full_path"
		exit 1
	fi
	awk 'BEGIN {max=0;read=0;band=1}
	{
		if(read==0 && $2 == "Band" && $3 == band){read=1}
		else if(read==1 && $2 == "Band" && $3 == band)
			{printf"Max_wfc_%d %.4f\n",band,max;read =0;band+=1;max=0}
		else if(read==1)
			{
				for(i=1;i<=NF;i++)
				{
					if(sqrt($i*$i)>max) {max=sqrt($i*$i)}
				}
			}
	}' "$full_path" >> "$props_result_file"
fi


#--------------------------------------------
# wave functions in LCAO basis
# echo "$has_lowf" # test out_wfc_lcao > 0
#--------------------------------------------
if ! test -z "$has_lowf" && [ "$has_lowf" == 1 ] && [ "${LIBRPA_PRODUCER_CONTRACT:-0}" != 1 ]; then
	if ! test -z "$gamma_only"  && [ $gamma_only == 1 ]; then
		wfc_cal=OUT.autotest/wf_nao.txt
		wfc_ref=wf_nao.txt.ref
	else  # multi-k point case
		if ! test -z "$out_app_flag"  && [ $out_app_flag == 0 ]; then
			wfc_name=wfk1g3_nao
			input_file=WFC/wfk1g3_nao
		else
			wfc_name=wfk2_nao
			input_file=wfk2_nao
		fi
	awk 'BEGIN {flag=999}
	{
		if($2 == "(band)") {flag=2;print $0}
		else if(flag>0) {flag-=1;print $0}
		else if(flag==0)
		{
			for(i=1;i<=NF/2;i++)
			{printf "%.10e ",sqrt( $(2*i)*$(2*i)+$(2*i-1)*$(2*i-1) )};
			printf "\n"
		}
		else {print $0}
	}' OUT.autotest/"$input_file".txt > OUT.autotest/"$wfc_name"_mod.txt
		wfc_cal=OUT.autotest/"$wfc_name"_mod.txt
		wfc_ref="$wfc_name"_mod.txt.ref
	fi

	python3 $COMPARE_SCRIPT $wfc_cal $wfc_ref 8 -abs 1
	echo "Compare_wfc_lcao_pass $?" >>$props_result_file
fi

}

run_cube_tail_props(){

#--------------------------------------------
# mulliken charge
#--------------------------------------------
if ! test -z "$out_mul"  && [ $out_mul == 1 ]; then
    python3 $COMPARE_SCRIPT mulliken.txt.ref OUT.autotest/mulliken.txt 3
	echo "Compare_mulliken_pass $?" >>$props_result_file
fi

#--------------------------------------------
# Process .cube files for:
# 1. get_wf/get_pchg calculation tag (LCAO)
# 2. out_wfc_norm/out_wfc_re_im/out_pchg (PW)
#--------------------------------------------
shopt -s nullglob
for pchg_ref in pchgi*.cube.ref; do
    pchg_cube=${pchg_ref%.ref}
    pchg_key=$(sanitize_result_key "${pchg_cube}_compare")
    record_compare_result "$props_result_file" "$pchg_key" "$pchg_ref" "OUT.autotest/$pchg_cube" 6
done
shopt -u nullglob

need_process_cube=false
# Check if this is a LCAO calculation with get_wf/get_pchg
if [ $calculation == "get_wf" ] || [ $calculation == "get_pchg" ]; then
    need_process_cube=true
fi
# Check if this is a PW calculation with out_wfc_norm/out_wfc_re_im
out_wfc_norm=$(get_input_key_value "out_wfc_norm" "INPUT")
out_wfc_re_im=$(get_input_key_value "out_wfc_re_im" "INPUT")
out_pchg=$(get_input_key_value "out_pchg" "INPUT")
if [ -n "$out_wfc_norm" ] || [ -n "$out_wfc_re_im" ] || [ -n "$out_pchg" ]; then
    need_process_cube=true
fi
# Process .cube files if needed
if [ "$need_process_cube" = true ]; then
    cubefiles=$(ls OUT.autotest/ | grep -E '.cube$')
    wavefunction_re_files=()

    if [ -z "$cubefiles" ]; then
        echo "Error: No .cube files found in OUT.autotest/"
        exit 1
    else
        for cube in $cubefiles; do
            if [[ "$cube" =~ ^wfi[0-9]+s[0-9]+(k[0-9]+)?re[.]cube$ ]]; then
                wavefunction_re_files+=("$cube")
                continue
            fi
            if [[ "$cube" =~ ^wfi[0-9]+s[0-9]+(k[0-9]+)?im[.]cube$ ]]; then
                continue
            fi
            total_chg=$(python3 "$CUBE_TOOL" integrate "OUT.autotest/$cube")
            echo "$cube $total_chg" >> $props_result_file
        done
    fi

    for cube in "${wavefunction_re_files[@]}"; do
        if [[ "$cube" =~ ^wfi([0-9]+)s([0-9]+)(k[0-9]+)?re[.]cube$ ]]; then
            band=${BASH_REMATCH[1]}
            spin=${BASH_REMATCH[2]}
            kpoint=${BASH_REMATCH[3]}
            state_prefix=${cube%re.cube}
            fingerprint_args=(
                "OUT.autotest/$cube"
                "OUT.autotest/${state_prefix}im.cube"
            )
            if [ "$nspin" = "4" ]; then
                if [ "$spin" != "1" ]; then
                    continue
                fi
                lower_prefix="wfi${band}s2${kpoint}"
                fingerprint_args+=(
                    "OUT.autotest/${lower_prefix}re.cube"
                    "OUT.autotest/${lower_prefix}im.cube"
                )
                result_prefix="wfi${band}${kpoint}_spinor_wfc_fp"
            else
                result_prefix="${state_prefix}_wfc_fp"
            fi

            if fingerprint=$(python3 "$CUBE_TOOL" fingerprint-wfc "${fingerprint_args[@]}"); then
                while read -r metric value; do
                    echo "${result_prefix}_${metric} $value" >> "$props_result_file"
                done <<< "$fingerprint"
            else
                echo "Error: Failed to generate wavefunction fingerprint for $state_prefix"
                exit 1
            fi
        fi
    done
fi

# Check the pointwise Pauli identities when all PW nspin=4 spinor outputs are available.
nspin=$(get_input_key_value "nspin" "INPUT")
if_separate_k=$(get_input_key_value "if_separate_k" "INPUT")
if [ "$nspin" = "4" ] && { [ "$if_separate_k" = "1" ] || [ "$if_separate_k" = "true" ]; } \
    && [ -n "$out_wfc_norm" ] && [ -n "$out_wfc_re_im" ] && [ -n "$out_pchg" ]; then
    if python3 "$CUBE_TOOL" check-spinor OUT.autotest; then
        echo "pw_spinor_cube_identity 0" >> "$props_result_file"
    else
        echo "pw_spinor_cube_identity 1" >> "$props_result_file"
    fi
fi

}
