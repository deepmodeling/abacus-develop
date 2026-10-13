#!/bin/bash

# mohan add 2025-05-03
# this compare script is used in different integrate tests
#
# Entry point: sources the shared library (tool paths, helper functions,
# one-time INPUT switch parsing) and then runs each property block below.
PROPS_SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$PROPS_SCRIPT_DIR/props_common.sh"
props_init "$1"
source "$PROPS_SCRIPT_DIR/props_basic.sh"
source "$PROPS_SCRIPT_DIR/props_mat.sh"
source "$PROPS_SCRIPT_DIR/props_cube.sh"
source "$PROPS_SCRIPT_DIR/props_ml.sh"
source "$PROPS_SCRIPT_DIR/props_tddft.sh"
source "$PROPS_SCRIPT_DIR/props_deepks.sh"

# Property collectors run in the same order as the original monolithic
# script so result files stay byte-identical; props_finalize() writes the
# trailing totaltimeref entry. Add a new category by dropping a
# props_<cat>.sh next to the others and calling run_<cat>_props here.
run_basic_props

#-------------------------------
# matrix collectors: out_dm1 (props_mat.sh)
#-------------------------------
run_mat_dm1_props

#-------------------------------
# cube collectors: out_pot/out_elf (props_cube.sh)
#-------------------------------
run_cube_pot_props

#-------------------------------
# matrix collectors (props_mat.sh): S/H(k)/XC/eband/H(R)/NPZ/r/T/dH
#-------------------------------
run_mat_props

#---------------------------------------
# cube collectors (props_cube.sh): chg/tau/LDOS/wfc real-space/PW/LCAO
#---------------------------------------
run_cube_props

#--------------------------------------------
# matrix collector: out_dm (props_mat.sh)
#--------------------------------------------
run_mat_dm_props

#--------------------------------------------
# cube collectors (props_cube.sh): mulliken/pchg/fingerprints/spinor
#--------------------------------------------
run_cube_tail_props

#--------------------------------------------
# ML collectors (props_ml.sh): MLKEDF descriptors
#--------------------------------------------
run_ml_descriptor_props

#--------------------------------------------
# basic collectors that run after ml: imp_sol
#--------------------------------------------
run_basic_props_post_ml

#--------------------------------------------
# ML collectors (props_ml.sh): RPA
#--------------------------------------------
run_ml_rpa_props

#--------------------------------------------
# DeePKS collectors (props_deepks.sh)
#--------------------------------------------
# Before the split this block ran as a separate `bash` process without -e,
# so a failing helper (get_sum_*.py, a missing deepks_desc.dat) left an
# empty value for the threshold check instead of aborting the collection.
# Keep that behavior with a subshell that disables errexit.
( set +e; run_deepks_props )

#--------------------------------------------
# basic collectors that run after deepks: symmetry
#--------------------------------------------
run_basic_props_post_deepks

#--------------------------------------------
# rt-TDDFT collectors (props_tddft.sh): current/efield/vecpot
#--------------------------------------------
run_tddft_props

#--------------------------------------------
# ML collectors (props_ml.sh): linear response excitations
#--------------------------------------------
run_ml_lr_props
#--------------------------------------------
# ML collectors (props_ml.sh): RDMFT energy terms
#--------------------------------------------
run_ml_rdmft_props

#--------------------------------------------
# basic collectors that run after rdmft: alllog
#--------------------------------------------
run_basic_props_post_rdmft

#--------------------------------------------
# trailing total-time entry
#--------------------------------------------
props_finalize "$1"
