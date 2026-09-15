#include "dfpt_metal.h"

#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

DFPT_Metal::DFPT_Metal()
{
}

DFPT_Metal::~DFPT_Metal()
{
}

void DFPT_Metal::init(double sigma, const std::string& smearing_type)
{
    ModuleBase::TITLE("DFPT_Metal", "init");
    ModuleBase::timer::start("DFPT_Metal", "init");
    sigma_ = sigma;
    smearing_type_ = smearing_type;
    ModuleBase::timer::end("DFPT_Metal", "init");
}

void DFPT_Metal::dfdeps(const ModuleBase::matrix& eig, double efermi, ModuleBase::matrix& dfdeps)
{
    ModuleBase::TITLE("DFPT_Metal", "dfdeps");
    ModuleBase::timer::start("DFPT_Metal", "dfdeps");
    // C4 interface reservation: the metallic DFPT branch (smearing
    // derivatives, Fermi-level shift dmu and the occupation-response part of
    // the first-order density) is intentionally NOT implemented in this
    // design-phase iteration; insulating systems only. The data members
    // (sigma_, smearing_type_) and the is_metal_/dmu_ slots of DFPT_PW_Data
    // are already in place for the future implementation.
    (void)eig;
    (void)efermi;
    (void)dfdeps;
    ModuleBase::WARNING_QUIT("DFPT_Metal", "metallic DFPT (dfdeps) is not supported in the design phase");
}

void DFPT_Metal::compute_dmu(int q_idx,
                             const psi::Psi<std::complex<double>>& psi,
                             const ModuleBase::matrix& wg,
                             const ModuleBase::matrix& dfdeps,
                             DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Metal", "compute_dmu");
    ModuleBase::timer::start("DFPT_Metal", "compute_dmu");
    (void)q_idx;
    (void)psi;
    (void)wg;
    (void)dfdeps;
    (void)data;
    ModuleBase::WARNING_QUIT("DFPT_Metal", "metallic DFPT (compute_dmu) is not supported in the design phase");
}

void DFPT_Metal::compute_drho_metal(int q_idx,
                                    const psi::Psi<std::complex<double>>& psi,
                                    const ModuleBase::matrix& wg,
                                    const ModuleBase::matrix& dfdeps,
                                    double dmu,
                                    DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Metal", "compute_drho_metal");
    ModuleBase::timer::start("DFPT_Metal", "compute_drho_metal");
    (void)q_idx;
    (void)psi;
    (void)wg;
    (void)dfdeps;
    (void)dmu;
    (void)data;
    ModuleBase::WARNING_QUIT("DFPT_Metal", "metallic DFPT (compute_drho_metal) is not supported in the design phase");
}

double DFPT_Metal::fd_dfdeps(double e, double efermi)
{
    ModuleBase::TITLE("DFPT_Metal", "fd_dfdeps");
    ModuleBase::timer::start("DFPT_Metal", "fd_dfdeps");
    (void)e;
    (void)efermi;
    ModuleBase::timer::end("DFPT_Metal", "fd_dfdeps");
    return 0.0;
}

double DFPT_Metal::gauss_dfdeps(double e, double efermi)
{
    ModuleBase::TITLE("DFPT_Metal", "gauss_dfdeps");
    ModuleBase::timer::start("DFPT_Metal", "gauss_dfdeps");
    (void)e;
    (void)efermi;
    ModuleBase::timer::end("DFPT_Metal", "gauss_dfdeps");
    return 0.0;
}

} // namespace ModuleDFPT
