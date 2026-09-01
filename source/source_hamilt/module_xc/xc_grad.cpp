// This file implements the public entry point of the gradient
// correction, XC_Functional::gradcorr, which drives a three-stage
// pipeline:
// 1. gradcorr_prepare_rho (xc_grad_prepare.cpp): builds the real- and
//      reciprocal-space density buffers and their gradients
// 2. gradcorr_xc_kernel (xc_grad_kernel.cpp): evaluates the GGA/meta-GGA
//      kernel on every grid point
// 3. gradcorr_assemble_vxc (xc_grad_assemble.cpp): assembles the
//      correction into etxc/vtxc/v/stress_gga
// The internal stage functions and their shared parameter/buffer
// structs are declared in xc_grad_internal.h.
// xc_grad_utils.cpp / xc_grad_wfc.cpp provide the auxiliary
// subroutines: grad_wfc, grad_rho, grad_dot, laplacian_rho and
// noncolin_rho.

#include "xc_functional.h"
#include "xc_grad_internal.h"
#include "source_base/timer.h"

void XC_Functional::gradcorr(
    double &etxc,
    double &vtxc,
    ModuleBase::matrix &v,
    const Charge* const chr,
    ModulePW::PW_Basis* rhopw,
    const UnitCell *ucell,
    std::vector<double> &stress_gga,
    const bool is_stress,
    const int nspin,
    const bool domag,
    const bool domag_z,
    const double hybrid_alpha_in,
    const double hse_omega_in)
{
    ModuleBase::TITLE("XC_Functional","gradcorr");

    if((func_type == 3 || func_type == 5) && nspin==4)
    {
        ModuleBase::WARNING_QUIT("gradcorr","meta-GGA has not been implemented for nspin = 4 yet");
    }

    if(func_type == 0 || func_type == 1)
    {
        return;
    }

    // No (semi-)local functional at all. `set_xc_type("HF")` sets func_type = 4 but leaves func_id empty,
    // Without this block, it will read back whatever the previous functional left in the cleared-but-not-freed
    // buffer (PBE, from set_xc_first_loop), producing a phantom XC energy and potential. (issue deepmodeling/abacus-develop#5404)
    if(func_id.empty())
    {
        // `Stress_Func::stress_gga` only guards on func_type (0/1) and unconditionally reads
        // stress_gga[0..8] afterwards, so hand back an explicit zero tensor rather than an
        // untouched (empty) vector.
        if(is_stress)
        {
            stress_gga.assign(9, 0.0);
        }
        return;
    }

    bool igcc_is_lyp = false;
    // func_id may hold a single entry (e.g. PBE0 -> {XC_HYB_GGA_XC_PBEH}), so guard the index.
    if( func_id.size() > 1 && func_id[1] == XC_GGA_C_LYP)
    {
        igcc_is_lyp = true;
    }

    int nspin0 = nspin;
    if(nspin==4)
    {
        nspin0 =1;
    }
    if(nspin==4&&(domag||domag_z))
    {
        nspin0 = 2;
    }

    assert(nspin0>0);
    const double fac = 1.0/ nspin0;

    // Use cached need_laplacian from set_xc_type
    bool need_laplacian = XC_Functional::get_need_laplacian();

    if(is_stress)
    {
        stress_gga.resize(9);
        for(int i=0;i<9;i++)
        {
            stress_gga[i] = 0.0;
        }
    }

    // Parameters and scratch buffers shared by the three pipeline stages.
    GradCorrParams params;
    params.chr = chr;
    params.rhopw = rhopw;
    params.ucell = ucell;
    params.nspin = nspin;
    params.nspin0 = nspin0;
    params.fac = fac;
    params.need_laplacian = need_laplacian;
    params.is_stress = is_stress;
    params.igcc_is_lyp = igcc_is_lyp;
    params.domag = domag;
    params.domag_z = domag_z;
    params.hybrid_alpha = hybrid_alpha_in;
    params.hse_omega = hse_omega_in;
    params.use_libxc = use_libxc;
    params.func_type = func_type;
    params.func_id = func_id;

    GradCorrBuffers buf;

    gradcorr_prepare_rho(params, v, buf);

    double vtxcgc = 0.0;
    double etxcgc = 0.0;

    gradcorr_xc_kernel(params, vtxcgc, etxcgc, stress_gga, v, buf);

    gradcorr_assemble_vxc(params, vtxcgc, etxcgc, vtxc, etxc, stress_gga, v, buf);

    return;
}
