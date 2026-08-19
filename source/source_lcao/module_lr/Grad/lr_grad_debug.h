#pragma once
// TEMPORARY DEBUG SCAFFOLDING -- ablation switches for the LR-KERNEL EXX entry points
// (the ones gated on `xc_kernel`, i.e. present in the `hf` tier but not in `rpa_at_hf`).
// These are gradient-only: none of them changes Omega, so the finite-difference reference
// stays fixed and the ablation is meaningful.
//   LRDBG_EXX_ZRK   Z-vector RHS, the K[D^X] term
//   LRDBG_EXX_KEXX  EDM, the W^X = 2K[D^X] term (`op_K_exx`)
//   LRDBG_EXX_FDMT  force, `cal_force_exx_dm_trans`
// Unset or empty => 1.0.
#include <cstdlib>
#include <string>
namespace LR_DBG
{
    inline double exx_scale(const char* key)
    {
        const std::string var = std::string("LRDBG_EXX_") + key;
        const char* v = std::getenv(var.c_str());
        return (v && *v) ? std::atof(v) : 1.0;
    }
}
