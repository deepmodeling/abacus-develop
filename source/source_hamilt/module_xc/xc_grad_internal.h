//==========================================================
// Internal declarations of the gradient-correction (gradcorr)
// pipeline. This header is included only by the xc_grad*.cpp
// implementation files; the public entry point
// XC_Functional::gradcorr stays in xc_functional.h.
//==========================================================
#ifndef XC_GRAD_INTERNAL_H
#define XC_GRAD_INTERNAL_H

#include "source_base/matrix.h"
#include "source_base/vector3.h"

#include <complex>
#include <vector>

class Charge;
class UnitCell;
namespace ModulePW
{
class PW_Basis;
}

/**
 * @brief Parameters shared by all stages of the gradcorr pipeline.
 *
 * Built once by XC_Functional::gradcorr and passed by const reference.
 * Packing the parameters into one struct keeps same-typed arguments
 * from being silently swapped and leaves stage signatures stable when
 * a new parameter is needed.
 */
struct GradCorrParams
{
    const Charge* chr;
    ModulePW::PW_Basis* rhopw;
    const UnitCell* ucell;
    int nspin;
    int nspin0; // number of spin channels of the (rho, |grad rho|) pair
    double fac; // 1.0 / nspin0
    bool need_laplacian;
    bool is_stress;
    bool igcc_is_lyp;
    bool domag;
    bool domag_z;
    double hybrid_alpha;
    double hse_omega;
    bool use_libxc;
    int func_type;            // 0:none, 1:lda, 2:gga, 3:mgga, 4:hybrid lda/gga, 5:hybrid mgga
    std::vector<int> func_id; // libxc id(s) of the functional
};

/**
 * @brief Scratch buffers shared by the gradcorr pipeline stages.
 *
 * All arrays are allocated by gradcorr_prepare_rho and consumed by
 * gradcorr_xc_kernel and gradcorr_assemble_vxc. The flat vsave/vgg
 * arrays are indexed as [is * nrxx + ir].
 */
struct GradCorrBuffers
{
    std::vector<double> rhotmp1;
    std::vector<double> rhotmp2;
    std::vector<std::complex<double>> rhogsum1;
    std::vector<std::complex<double>> rhogsum2;
    std::vector<ModuleBase::Vector3<double>> gdr1;
    std::vector<ModuleBase::Vector3<double>> gdr2;
    std::vector<ModuleBase::Vector3<double>> h1;
    std::vector<ModuleBase::Vector3<double>> h2;
    std::vector<double> neg;
    std::vector<double> vsave;
    std::vector<double> vgg;
    std::vector<double> lapl1;
    std::vector<double> lapl2;
    std::vector<double> vlapl_arr1;
    std::vector<double> vlapl_arr2;
};

/**
 * @brief First stage: build the real- and reciprocal-space density
 *        buffers and their gradients for every spin channel.
 */
void gradcorr_prepare_rho(const GradCorrParams& params,
                          ModuleBase::matrix& v,
                          GradCorrBuffers& buf);

/**
 * @brief Second stage: evaluate the GGA/meta-GGA kernel on every grid
 *        point and accumulate vtxcgc, etxcgc, stress_gga and v.
 */
void gradcorr_xc_kernel(const GradCorrParams& params,
                        double& vtxcgc,
                        double& etxcgc,
                        std::vector<double>& stress_gga,
                        ModuleBase::matrix& v,
                        GradCorrBuffers& buf);

/**
 * @brief Final stage: assemble the gradient-correction contribution
 *        into vtxc/etxc/stress_gga/v and restore the nspin = 4 case.
 */
void gradcorr_assemble_vxc(const GradCorrParams& params,
                           double vtxcgc,
                           double etxcgc,
                           double& vtxc,
                           double& etxc,
                           std::vector<double>& stress_gga,
                           ModuleBase::matrix& v,
                           GradCorrBuffers& buf);

#endif // XC_GRAD_INTERNAL_H
