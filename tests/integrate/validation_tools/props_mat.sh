#!/bin/bash

# Matrix/operator property collectors (S, H(k), H(R), XC, r/T/dH(R), DM, ...).
#
# Sourced by catch_properties.sh after props_common.sh; relies on the global
# switch variables set by props_init() and appends result lines to
# $props_result_file.
#
# The original monolithic script interleaved these checks with the cube
# collectors (props_cube.sh). To preserve the exact emission order of result
# lines this module exposes three hooks, called by the entry point at the
# matching positions:
#   run_mat_dm1_props() - out_dm1 block (before the pot/ELF cube blocks)
#   run_mat_props()     - everything from S(k) through dH(k) term matrices
#   run_mat_dm_props()  - out_dm block (between wfc_lcao and mulliken)

# Compare an H(R) CSR file with its reference.
# The "# spin index, E_Fermi = <x> eV" header annotation is printed with
# setprecision(6), so the matrix body is compared at 1e-8 with the annotation
# stripped, and the E_Fermi values are compared separately at 1e-3 eV.
# A missing file is reported as a failure (process substitution would
# otherwise hand CompareFile.py an empty stream and pass silently).
compare_hr_csr(){
	local ref=$1
	local cal=$2
	local strip='s/, E_Fermi = [0-9.eE+-]* eV//'
	local pick='s/.*E_Fermi = \([0-9.eE+-]*\) eV.*/\1/p'
	if ! test -f "$ref" || ! test -f "$cal"; then
		echo "Error: can not find file $ref or $cal"
		return 1
	fi
	python3 "$COMPARE_SCRIPT" <(sed "$strip" "$ref") <(sed "$strip" "$cal") 8 || return 1
	if [ "$(sed -n "$pick" "$ref" | wc -l)" -ne "$(sed -n "$pick" "$cal" | wc -l)" ]; then
		echo "Error: E_Fermi annotation count differs between $ref and $cal"
		return 1
	fi
	python3 "$COMPARE_SCRIPT" <(sed -n "$pick" "$ref") <(sed -n "$pick" "$cal") 3
}

run_mat_dm1_props(){

#-------------------------------
# echo $out_dm1
#-------------------------------
if ! test -z "$out_dm1"  && [  $out_dm1 == 1 ]; then
	dm1ref=dmr_nao.csr.ref
	dm1cal=OUT.autotest/dmr_nao.csr
	record_compare_result "$props_result_file" "CompareDM1_pass" "$dm1ref" "$dm1cal" 8
fi

}

run_mat_props(){

#-------------------------------
# Overlap matrix
# calculation == get_s
#-------------------------------
if [ "$calculation" == "get_s" ]; then
	sref=sr_nao.csr.ref
	scal=OUT.autotest/sr_nao.csr
	record_compare_result "$props_result_file" "CompareS_pass" "$sref" "$scal" 8
fi

#-------------------------------
# Partial band structure
# echo $out_pband
#-------------------------------
if ! test -z "$out_pband"  && [  $out_pband == 1 ]; then
	orbref=refOrbital
	orbcal=OUT.autotest/Orbital
	record_compare_result "$props_result_file" "CompareOrb_pass" "$orbref" "$orbcal" 8
fi

#-------------------------------
# Wannier90 information
# echo $toW90
#-------------------------------
if ! test -z "$toW90"  && [  $toW90 == 1 ]; then
	amnref=diamond.amn
	amncal=OUT.autotest/diamond.amn
	mmnref=diamond.mmn
	mmncal=OUT.autotest/diamond.mmn
	eigref=diamond.eig
	eigcal=OUT.autotest/diamond.eig
	sed -i '1d' "$amncal"
	sed -i '1d' "$mmncal"
	record_compare_result "$props_result_file" "CompareAMN_pass" "$amnref" "$amncal" 1 8
	record_compare_result "$props_result_file" "CompareMMN_pass" "$mmnref" "$mmncal" 1 8
	record_compare_result "$props_result_file" "CompareEIG_pass" "$eigref" "$eigcal" 8
fi

#-------------------------------
# Total DOS
# echo total_dos
# echo $has_band
#-------------------------------
if ! test -z "$has_band"  && [  $has_band == 1 ]; then
	bandref=band.txt.ref
	bandcal=OUT.autotest/band.txt
	record_compare_result "$props_result_file" "CompareBand_pass" "$bandref" "$bandcal" 8
fi

#--------------------------------
# Hamiltonian and overlap matrix
# echo $has_hs
#--------------------------------
if ! test -z "$has_hs"  && [ $has_hs == 1 ]; then
    if ! test -z "$gamma_only"  && [ $gamma_only == 1 ]; then
        # ========== Γ-point (single k-point) calculation ==========
        if ! test -z "$nspin" && [ $nspin == 2 ]; then
            # nspin=2 (spin-polarized): compare hks1 + hks2 Hamiltonian + sk overlap matrix
            h1ref=hks1_nao.txt.ref
            h1cal=OUT.autotest/hks1_nao.txt
            h2ref=hks2_nao.txt.ref
            h2cal=OUT.autotest/hks2_nao.txt
            sref=sk_nao.txt.ref
            scal=OUT.autotest/sk_nao.txt
            # Compare Hamiltonian matrix for spin 1
            record_compare_result "$props_result_file" "CompareH1_pass" "$h1ref" "$h1cal" 6
            # Compare Hamiltonian matrix for spin 2
            record_compare_result "$props_result_file" "CompareH2_pass" "$h2ref" "$h2cal" 6
            # Compare overlap matrix
            record_compare_result "$props_result_file" "CompareS_pass" "$sref" "$scal" 8
        elif ! test -z "$nspin" && [ $nspin == 4 ]; then
            # nspin=4 : do nothing, only matching condition without any operation
            true
        else
            # nspin=1 (non-spin-polarized, default case): compare single hk + sk matrix set
            href=hk_nao.txt.ref
            hcal=OUT.autotest/hk_nao.txt
            sref=sk_nao.txt.ref
            scal=OUT.autotest/sk_nao.txt
            # Compare Hamiltonian matrix
            record_compare_result "$props_result_file" "CompareH_pass" "$href" "$hcal" 6
            # Compare overlap matrix
            record_compare_result "$props_result_file" "CompareS_pass" "$sref" "$scal" 8
        fi
    else
        # ========== Multiple k-points calculation ==========
        if ! test -z "$nspin" && [ $nspin == 2 ]; then
            # nspin=2 (spin-polarized): compare spin-up/spin-down H(k) and S(k) at the second k-point
            h1ref=hk2s1_nao.txt.ref
            h1cal=OUT.autotest/hk2s1_nao.txt
            h2ref=hk2s2_nao.txt.ref
            h2cal=OUT.autotest/hk2s2_nao.txt
            sref=sk2_nao.txt.ref
            scal=OUT.autotest/sk2_nao.txt
            # Compare Hamiltonian matrix for spin 1
            record_compare_result "$props_result_file" "CompareH1_pass" "$h1ref" "$h1cal" 6
            # Compare Hamiltonian matrix for spin 2
            record_compare_result "$props_result_file" "CompareH2_pass" "$h2ref" "$h2cal" 6
            # Compare overlap matrix
            record_compare_result "$props_result_file" "CompareS_pass" "$sref" "$scal" 8
        elif ! test -z "$nspin" && [ $nspin == 4 ]; then
            # nspin=4 : do nothing, only matching condition without any operation
            true
        else
            # nspin=1 (non-spin-polarized, default case): compare single hk2 + sk2 matrix set
            href=hk2_nao.txt.ref
            hcal=OUT.autotest/hk2_nao.txt
            sref=sk2_nao.txt.ref
            scal=OUT.autotest/sk2_nao.txt
            # Compare Hamiltonian matrix
            record_compare_result "$props_result_file" "CompareH_pass" "$href" "$hcal" 6
            # Compare overlap matrix
            record_compare_result "$props_result_file" "CompareS_pass" "$sref" "$scal" 8
        fi
    fi
elif ! test -z "$has_hs" && [ $has_hs == 2 ]; then
    HSK_BINARY_COMPARE="$PROPS_TOOLS_DIR/compare_hsk_binary.py"
    if ! test -z "$gamma_only" && [ $gamma_only == 1 ]; then
        HSK_TEXT_REFERENCE_DIR="../scf_out_hk"
        run_and_record_status "$props_result_file" "CompareH_pass" python3 "$HSK_BINARY_COMPARE" \
            OUT.autotest/hk_nao.dat "$HSK_TEXT_REFERENCE_DIR/hk_nao.txt.ref" real 3
        run_and_record_status "$props_result_file" "CompareS_pass" python3 "$HSK_BINARY_COMPARE" \
            OUT.autotest/sk_nao.dat "$HSK_TEXT_REFERENCE_DIR/sk_nao.txt.ref" real 3
    else
        HSK_TEXT_REFERENCE_DIR="../scf_out_hsk"
        run_and_record_status "$props_result_file" "CompareH_pass" python3 "$HSK_BINARY_COMPARE" \
            OUT.autotest/hk2_nao.dat "$HSK_TEXT_REFERENCE_DIR/hk2_nao.txt.ref" complex 3
        run_and_record_status "$props_result_file" "CompareS_pass" python3 "$HSK_BINARY_COMPARE" \
            OUT.autotest/sk2_nao.dat "$HSK_TEXT_REFERENCE_DIR/sk2_nao.txt.ref" complex 3
    fi
fi

#--------------------------------
# exchange-correlation potential
#--------------------------------
if ! test -z "$has_xc"  && [  $has_xc == 1 ]; then
	if ! test -z "$gamma_only"  && [ $gamma_only == 1 ]; then
			xcref=vxc_nao.txt.ref
			xccal=OUT.autotest/vxc_nao.txt
	else
			xcref=vxck2_nao.txt.ref
			xccal=OUT.autotest/vxck2_nao.txt
	fi
	oeref=vxc_out.ref
	oecal=OUT.autotest/vxc_out.dat
	record_compare_result "$props_result_file" "CompareVXC_pass" "$xcref" "$xccal" 4
	record_compare_result "$props_result_file" "CompareOrbXC_pass" "$oeref" "$oecal" 5
fi

#--------------------------------
# exchange-correlation potential
#--------------------------------
if ! test -z "$has_xc2"  && [  $has_xc2 == 1 ]; then
	xc2ref=Vxc_R_spin0.ref
	xc2cal=OUT.autotest/Vxc_R_spin0.csr
	record_compare_result "$props_result_file" "CompareVXC_R_pass" "$xc2ref" "$xc2cal" 8
fi

#--------------------------------
# separate terms in band enegy
#--------------------------------
if ! test -z "$has_eband_separate"  && [  $has_eband_separate == 1 ]; then
	ekref=kinetic_out.ref
	ekcal=OUT.autotest/kinetic_out.dat
	record_compare_result "$props_result_file" "CompareOrbKinetic_pass" "$ekref" "$ekcal" 4
	vlref=vpp_local_out.ref
	vlcal=OUT.autotest/vpp_local_out.dat
	record_compare_result "$props_result_file" "CompareOrbVL_pass" "$vlref" "$vlcal" 4
	vnlref=vpp_nonlocal_out.ref
	vnlcal=OUT.autotest/vpp_nonlocal_out.dat
	record_compare_result "$props_result_file" "CompareOrbVNL_pass" "$vnlref" "$vnlcal" 4
	vhref=vhartree_out.ref
	vhcal=OUT.autotest/vhartree_out.dat
	record_compare_result "$props_result_file" "CompareOrbVHartree_pass" "$vhref" "$vhcal" 4
fi

#-----------------------------------
# Hamiltonian and overlap matrices
#-----------------------------------
#echo $has_hs2
if ! test -z "$has_hs2"  && [  $has_hs2 == 1 ]; then
    run_and_record_status "$props_result_file" "CompareHR_pass" \
        compare_hr_csr hrs1_nao.csr.ref OUT.autotest/hrs1_nao.csr
    if ! test -z "$nspin" && [ "$nspin" -eq 2 ]; then
        run_and_record_status "$props_result_file" "CompareHR2_pass" \
            compare_hr_csr hrs2_nao.csr.ref OUT.autotest/hrs2_nao.csr
    fi
    record_compare_result "$props_result_file" "CompareSR_pass" \
        "sr_nao.csr.ref" "OUT.autotest/sr_nao.csr" 8
elif ! test -z "$has_hs2" && [ "$has_hs2" == 2 ]; then
    HSR_BINARY_COMPARE="$PROPS_TOOLS_DIR/compare_hsr_binary.py"
    run_and_record_status "$props_result_file" "CompareHR_pass" python3 "$HSR_BINARY_COMPARE" \
        OUT.autotest/hrs1_nao.dat hrs1_nao.csr.ref real 4
    if ! test -z "$nspin" && [ "$nspin" -eq 2 ]; then
        run_and_record_status "$props_result_file" "CompareHR2_pass" python3 "$HSR_BINARY_COMPARE" \
            OUT.autotest/hrs2_nao.dat hrs2_nao.csr.ref real 4
    fi
    run_and_record_status "$props_result_file" "CompareSR_pass" python3 "$HSR_BINARY_COMPARE" \
        OUT.autotest/sr_nao.dat sr_nao.csr.ref real 4
fi

#-----------------------------------
# H(R), S(R), and DM(R) matrices in NPZ format
#-----------------------------------
if { ! test -z "$out_hsr" && [ "$out_hsr" == 3 ]; } || { ! test -z "$out_hsr_npz" && [ "$out_hsr_npz" == 1 ]; }; then
    run_and_record_status "$props_result_file" "OutputSRNPZ_pass" test -f OUT.autotest/sr_nao.npz
fi

if { ! test -z "$out_hr_npz" && [ "$out_hr_npz" == 1 ]; } \
    || { ! test -z "$out_hsr" && [ "$out_hsr" == 3 ]; } \
    || { ! test -z "$out_hsr_npz" && [ "$out_hsr_npz" == 1 ]; }; then
    run_and_record_status "$props_result_file" "OutputHRNPZ_pass" test -f OUT.autotest/hrs1_nao.npz
fi

if ! test -z "$out_dm_npz" && [ "$out_dm_npz" == 1 ]; then
    run_and_record_status "$props_result_file" "OutputDMNPZ_pass" test -f OUT.autotest/output_DM0.npz
fi

#-----------------------------------
#  <psi_i0 | r | psi_jR> matrix
#-----------------------------------
#echo $has_mat_r
if ! test -z "$has_mat_r"  && [  $has_mat_r == 1 ]; then
    record_compare_result "$props_result_file" "ComparerR_pass" \
        "rr_nao.txt.ref" "OUT.autotest/rr_nao.txt" 8
fi

#-----------------------------------
#  <psi_i0 | T | psi_jR> matrix
#-----------------------------------
#echo $has_mat_t
if ! test -z "$has_mat_t"  && [  $has_mat_t == 1 ]; then
    record_compare_result "$props_result_file" "ComparerTR_pass" \
        "tr_nao.csr.ref" "OUT.autotest/tr_nao.csr" 8
fi

#-----------------------------------
#  Asynchronous overlap matrix for Hefei-NAMD
#-----------------------------------
#echo $has_mat_syns
if ! test -z "$has_mat_syns"  && [  $has_mat_syns == 1 ]; then
    record_compare_result "$props_result_file" "CompareSYNS_pass" \
        "syns_nao.csr.ref" "OUT.autotest/syns_nao.csr" 8
fi

#-----------------------------------
#  <psi_i0 | H | dpsi_jR> matrix
#-----------------------------------
#echo $has_mat_dh
if ! test -z "$has_mat_dh"  && [  $has_mat_dh == 1 ] && [ $gamma_only != 1 ]; then
    record_compare_result "$props_result_file" "ComparerdHRx_pass" \
        "dhrxs1_nao.csr.ref" "OUT.autotest/dhrxs1_nao.csr" 8
    record_compare_result "$props_result_file" "ComparerdHRy_pass" \
        "dhrys1_nao.csr.ref" "OUT.autotest/dhrys1_nao.csr" 8
    record_compare_result "$props_result_file" "ComparerdHRz_pass" \
        "dhrzs1_nao.csr.ref" "OUT.autotest/dhrzs1_nao.csr" 8
fi

#-----------------------------------
#  d <psi_i0 | H | psi_j>(k) matrix
#-----------------------------------
#echo $has_mat_dh_terms
if ! test -z "$has_mat_dh"  && [  $has_mat_dh == 1 ]; then
    shopt -s nullglob
    for reffile in dhk_ref/*.txt; do
        fname=$(basename "$reffile")
        key=$(sanitize_result_key "Compare_${fname%.txt}")
        record_compare_result "$props_result_file" "${key}_pass" "$reffile" "OUT.autotest/$fname" 8
    done
    shopt -u nullglob
fi

}

run_mat_dm_props(){

#--------------------------------------------
# density matrix information
#--------------------------------------------
if ! test -z "$out_dm"  && [ $out_dm == 1 ]; then
	dmfile=OUT.autotest/dm_nao.txt
	dmref=dm_nao.txt.ref
	record_compare_result "$props_result_file" "DM_different" "$dmref" "$dmfile" 5
fi

}
