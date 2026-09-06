#include "dfpt_pw_data.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include <cmath>

namespace ModuleDFPT
{

DFPT_PW_Data::DFPT_PW_Data()
{
}

DFPT_PW_Data::~DFPT_PW_Data()
{
    clean();
}

void DFPT_PW_Data::init(ModuleCell::QList* qlist,
                        int nk,
                        int nbands,
                        int npw_max,
                        int nrxx,
                        int nspin,
                        int nat,
                        const Plus_U_Base* dftu)
{
    ModuleBase::TITLE("DFPT_PW_Data", "init");
    ModuleBase::timer::start("DFPT_PW_Data", "init");
    qlist_ = qlist;
    nk_ = nk;
    nbands_ = nbands;
    npw_max_ = npw_max;
    nrxx_ = nrxx;
    nspin_ = nspin;
    nat_ = nat;
    dftu_ = dftu;

    allocate_memory();
    is_initialized_ = true;
    ModuleBase::timer::end("DFPT_PW_Data", "init");
}

void DFPT_PW_Data::clean()
{
    ModuleBase::TITLE("DFPT_PW_Data", "clean");
    ModuleBase::timer::start("DFPT_PW_Data", "clean");
    deallocate_memory();
    is_initialized_ = false;
    ModuleBase::timer::end("DFPT_PW_Data", "clean");
}

bool DFPT_PW_Data::u_active() const
{
    ModuleBase::TITLE("DFPT_PW_Data", "u_active");
    ModuleBase::timer::start("DFPT_PW_Data", "u_active");
    // a usable provider has its occupation matrices initialized (the ground
    // state does this when DFT+U actually runs); a wired provider without
    // them (e.g. a default-constructed reservation) stays inactive.
    ModuleBase::timer::end("DFPT_PW_Data", "u_active");
    return with_u() && dftu_->is_occ_mat_initialized();
}

void DFPT_PW_Data::set_docc(int q_idx, const std::vector<std::complex<double>>& occ)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_docc");
    ModuleBase::timer::start("DFPT_PW_Data", "set_docc");
    if (q_idx < 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_docc");
        return;
    }
    if (q_idx >= static_cast<int>(docc_.size()))
    {
        docc_.resize(q_idx + 1);
    }
    docc_[q_idx] = occ;
    ModuleBase::timer::end("DFPT_PW_Data", "set_docc");
}

void DFPT_PW_Data::set_vsc_r(int atom_idx, int dir, const std::vector<std::complex<double>>& v)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_vsc_r");
    ModuleBase::timer::start("DFPT_PW_Data", "set_vsc_r");
    if (atom_idx < 0 || dir < 0 || dir >= 3)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_vsc_r");
        return;
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot >= vsc_r_.size())
    {
        vsc_r_.resize(slot + 1);
    }
    vsc_r_[slot] = v;
    ModuleBase::timer::end("DFPT_PW_Data", "set_vsc_r");
}

std::vector<std::complex<double>> DFPT_PW_Data::get_vsc_r(int atom_idx, int dir) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_vsc_r");
    ModuleBase::timer::start("DFPT_PW_Data", "get_vsc_r");
    if (atom_idx < 0 || dir < 0 || dir >= 3)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_vsc_r");
        return std::vector<std::complex<double>>();
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot < vsc_r_.size())
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_vsc_r");
        return vsc_r_[slot];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_vsc_r");
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_dpsi_disp(int atom_idx,
                                 int dir,
                                 const std::vector<std::vector<std::vector<std::complex<double>>>>& d)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dpsi_disp");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dpsi_disp");
    if (atom_idx < 0 || dir < 0 || dir >= 3)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_dpsi_disp");
        return;
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot >= dpsi_disp_.size())
    {
        dpsi_disp_.resize(slot + 1);
    }
    dpsi_disp_[slot] = d;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dpsi_disp");
}

std::vector<std::vector<std::vector<std::complex<double>>>> DFPT_PW_Data::get_dpsi_disp(int atom_idx, int dir) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dpsi_disp");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dpsi_disp");
    if (atom_idx < 0 || dir < 0 || dir >= 3)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dpsi_disp");
        return std::vector<std::vector<std::vector<std::complex<double>>>>();
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot < dpsi_disp_.size())
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dpsi_disp");
        return dpsi_disp_[slot];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_dpsi_disp");
    return std::vector<std::vector<std::vector<std::complex<double>>>>();
}

void DFPT_PW_Data::set_pos_resp(int dir, const std::vector<std::vector<std::vector<std::complex<double>>>>& y)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_pos_resp");
    ModuleBase::timer::start("DFPT_PW_Data", "set_pos_resp");
    if (dir < 0 || dir >= 3)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_pos_resp");
        return;
    }
    if (pos_resp_.size() < 3)
    {
        pos_resp_.resize(3);
    }
    pos_resp_[dir] = y;
    ModuleBase::timer::end("DFPT_PW_Data", "set_pos_resp");
}

std::vector<std::vector<std::vector<std::complex<double>>>> DFPT_PW_Data::get_pos_resp(int dir) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_pos_resp");
    ModuleBase::timer::start("DFPT_PW_Data", "get_pos_resp");
    if (dir < 0 || dir >= 3 || dir >= static_cast<int>(pos_resp_.size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_pos_resp");
        return std::vector<std::vector<std::vector<std::complex<double>>>>();
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_pos_resp");
    return pos_resp_[dir];
}

void DFPT_PW_Data::set_dpsi_efield(int dir, const std::vector<std::vector<std::vector<std::complex<double>>>>& d)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dpsi_efield");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dpsi_efield");
    if (dir < 0 || dir >= 3)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_dpsi_efield");
        return;
    }
    if (dpsi_efield_.size() < 3)
    {
        dpsi_efield_.resize(3);
    }
    dpsi_efield_[dir] = d;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dpsi_efield");
}

std::vector<std::vector<std::vector<std::complex<double>>>> DFPT_PW_Data::get_dpsi_efield(int dir) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dpsi_efield");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dpsi_efield");
    if (dir < 0 || dir >= 3 || dir >= static_cast<int>(dpsi_efield_.size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dpsi_efield");
        return std::vector<std::vector<std::vector<std::complex<double>>>>();
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_dpsi_efield");
    return dpsi_efield_[dir];
}

std::vector<std::complex<double>> DFPT_PW_Data::get_docc(int q_idx) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_docc");
    ModuleBase::timer::start("DFPT_PW_Data", "get_docc");
    if (q_idx >= 0 && q_idx < static_cast<int>(docc_.size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_docc");
        return docc_[q_idx];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_docc");
    return std::vector<std::complex<double>>();
}

int DFPT_PW_Data::get_nq() const
{
    return qlist_->get_nq();
}

ModuleBase::Vector3<double> DFPT_PW_Data::get_qvec(int q_idx) const
{
    return qlist_->get_q(q_idx);
}

int DFPT_PW_Data::get_nirr(int q_idx) const
{
    return qlist_->get_nirr(q_idx);
}

std::vector<int> DFPT_PW_Data::get_irrep_modes(int q_idx, int irrep) const
{
    return qlist_->get_irrep_modes(q_idx, irrep);
}

void DFPT_PW_Data::set_dpsi(int q_idx, int k_idx, int band_idx, const std::vector<std::complex<double>>& psi)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dpsi");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dpsi");
    if (q_idx < 0 || k_idx < 0 || band_idx < 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_dpsi");
        return;
    }
    if (q_idx >= static_cast<int>(dpsi_.size()))
    {
        dpsi_.resize(q_idx + 1);
    }
    if (k_idx >= static_cast<int>(dpsi_[q_idx].size()))
    {
        dpsi_[q_idx].resize(k_idx + 1);
    }
    if (band_idx >= static_cast<int>(dpsi_[q_idx][k_idx].size()))
    {
        dpsi_[q_idx][k_idx].resize(band_idx + 1);
    }
    dpsi_[q_idx][k_idx][band_idx] = psi;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dpsi");
}

std::vector<std::complex<double>> DFPT_PW_Data::get_dpsi(int q_idx, int k_idx, int band_idx) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dpsi");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dpsi");
    if (q_idx >= 0 && k_idx >= 0 && band_idx >= 0 && q_idx < static_cast<int>(dpsi_.size())
        && k_idx < static_cast<int>(dpsi_[q_idx].size()) && band_idx < static_cast<int>(dpsi_[q_idx][k_idx].size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dpsi");
        return dpsi_[q_idx][k_idx][band_idx];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_dpsi");
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_converged(int q_idx, int irrep, bool flag)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_converged");
    ModuleBase::timer::start("DFPT_PW_Data", "set_converged");
    converged_[std::make_pair(q_idx, irrep)] = flag;
    ModuleBase::timer::end("DFPT_PW_Data", "set_converged");
}

bool DFPT_PW_Data::get_converged(int q_idx, int irrep) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_converged");
    ModuleBase::timer::start("DFPT_PW_Data", "get_converged");
    const auto it = converged_.find(std::make_pair(q_idx, irrep));
    ModuleBase::timer::end("DFPT_PW_Data", "get_converged");
    return it != converged_.end() ? it->second : false;
}

void DFPT_PW_Data::add_residual(int q_idx, int irrep, double r)
{
    ModuleBase::TITLE("DFPT_PW_Data", "add_residual");
    ModuleBase::timer::start("DFPT_PW_Data", "add_residual");
    residuals_[std::make_pair(q_idx, irrep)].push_back(r);
    ModuleBase::timer::end("DFPT_PW_Data", "add_residual");
}

std::vector<double> DFPT_PW_Data::get_residuals(int q_idx, int irrep) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_residuals");
    ModuleBase::timer::start("DFPT_PW_Data", "get_residuals");
    const auto it = residuals_.find(std::make_pair(q_idx, irrep));
    ModuleBase::timer::end("DFPT_PW_Data", "get_residuals");
    return it != residuals_.end() ? it->second : std::vector<double>();
}

void DFPT_PW_Data::set_current_iter(int q_idx, int irrep, int iter)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_current_iter");
    ModuleBase::timer::start("DFPT_PW_Data", "set_current_iter");
    current_iter_[std::make_pair(q_idx, irrep)] = iter;
    ModuleBase::timer::end("DFPT_PW_Data", "set_current_iter");
}

int DFPT_PW_Data::get_current_iter(int q_idx, int irrep) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_current_iter");
    ModuleBase::timer::start("DFPT_PW_Data", "get_current_iter");
    const auto it = current_iter_.find(std::make_pair(q_idx, irrep));
    ModuleBase::timer::end("DFPT_PW_Data", "get_current_iter");
    return it != current_iter_.end() ? it->second : 0;
}

void DFPT_PW_Data::set_drho_r(int q_idx, int spin, const std::vector<double>& rho)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_drho_r");
    ModuleBase::timer::start("DFPT_PW_Data", "set_drho_r");
    if (q_idx < 0 || spin < 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_drho_r");
        return;
    }
    if (q_idx >= static_cast<int>(drho_r_.size()))
    {
        drho_r_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(drho_r_[q_idx].size()))
    {
        drho_r_[q_idx].resize(spin + 1);
    }
    drho_r_[q_idx][spin] = rho;
    ModuleBase::timer::end("DFPT_PW_Data", "set_drho_r");
}

std::vector<double> DFPT_PW_Data::get_drho_r(int q_idx, int spin) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_drho_r");
    ModuleBase::timer::start("DFPT_PW_Data", "get_drho_r");
    if (q_idx >= 0 && spin >= 0 && q_idx < static_cast<int>(drho_r_.size())
        && spin < static_cast<int>(drho_r_[q_idx].size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_drho_r");
        return drho_r_[q_idx][spin];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_drho_r");
    return std::vector<double>();
}

void DFPT_PW_Data::set_drho_g(int q_idx, int spin, const std::vector<std::complex<double>>& rho)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_drho_g");
    ModuleBase::timer::start("DFPT_PW_Data", "set_drho_g");
    if (q_idx < 0 || spin < 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_drho_g");
        return;
    }
    if (q_idx >= static_cast<int>(drho_g_.size()))
    {
        drho_g_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(drho_g_[q_idx].size()))
    {
        drho_g_[q_idx].resize(spin + 1);
    }
    drho_g_[q_idx][spin] = rho;
    ModuleBase::timer::end("DFPT_PW_Data", "set_drho_g");
}

std::vector<std::complex<double>> DFPT_PW_Data::get_drho_g(int q_idx, int spin) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_drho_g");
    ModuleBase::timer::start("DFPT_PW_Data", "get_drho_g");
    if (q_idx >= 0 && spin >= 0 && q_idx < static_cast<int>(drho_g_.size())
        && spin < static_cast<int>(drho_g_[q_idx].size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_drho_g");
        return drho_g_[q_idx][spin];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_drho_g");
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_dv_r(int q_idx, int spin, const std::vector<double>& v)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dv_r");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dv_r");
    if (q_idx < 0 || spin < 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_dv_r");
        return;
    }
    if (q_idx >= static_cast<int>(dv_r_.size()))
    {
        dv_r_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(dv_r_[q_idx].size()))
    {
        dv_r_[q_idx].resize(spin + 1);
    }
    dv_r_[q_idx][spin] = v;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dv_r");
}

std::vector<double> DFPT_PW_Data::get_dv_r(int q_idx, int spin) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dv_r");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dv_r");
    if (q_idx >= 0 && spin >= 0 && q_idx < static_cast<int>(dv_r_.size())
        && spin < static_cast<int>(dv_r_[q_idx].size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dv_r");
        return dv_r_[q_idx][spin];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_dv_r");
    return std::vector<double>();
}

void DFPT_PW_Data::set_dv_recip_c(int q_idx, int spin, const std::vector<std::complex<double>>& v)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dv_recip_c");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dv_recip_c");
    if (q_idx < 0 || spin < 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_dv_recip_c");
        return;
    }
    if (q_idx >= static_cast<int>(dv_recip_c_.size()))
    {
        dv_recip_c_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(dv_recip_c_[q_idx].size()))
    {
        dv_recip_c_[q_idx].resize(spin + 1);
    }
    dv_recip_c_[q_idx][spin] = v;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dv_recip_c");
}

std::vector<std::complex<double>> DFPT_PW_Data::get_dv_recip_c(int q_idx, int spin) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dv_recip_c");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dv_recip_c");
    if (q_idx >= 0 && spin >= 0 && q_idx < static_cast<int>(dv_recip_c_.size())
        && spin < static_cast<int>(dv_recip_c_[q_idx].size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dv_recip_c");
        return dv_recip_c_[q_idx][spin];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_dv_recip_c");
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_dv_rc(int q_idx, int spin, const std::vector<std::complex<double>>& v)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dv_rc");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dv_rc");
    if (q_idx < 0 || spin < 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_dv_rc");
        return;
    }
    if (q_idx >= static_cast<int>(dv_rc_.size()))
    {
        dv_rc_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(dv_rc_[q_idx].size()))
    {
        dv_rc_[q_idx].resize(spin + 1);
    }
    dv_rc_[q_idx][spin] = v;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dv_rc");
}

std::vector<std::complex<double>> DFPT_PW_Data::get_dv_rc(int q_idx, int spin) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dv_rc");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dv_rc");
    if (q_idx >= 0 && spin >= 0 && q_idx < static_cast<int>(dv_rc_.size())
        && spin < static_cast<int>(dv_rc_[q_idx].size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dv_rc");
        return dv_rc_[q_idx][spin];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_dv_rc");
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_dynmat(int q_idx, const ModuleBase::ComplexMatrix& dm)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dynmat");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dynmat");
    if (q_idx >= static_cast<int>(dynmat_.size()))
    {
        dynmat_.resize(q_idx + 1);
    }
    dynmat_[q_idx] = dm;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dynmat");
}

ModuleBase::ComplexMatrix DFPT_PW_Data::get_dynmat(int q_idx) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dynmat");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dynmat");
    if (q_idx < static_cast<int>(dynmat_.size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dynmat");
        return dynmat_[q_idx];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_dynmat");
    return ModuleBase::ComplexMatrix();
}

void DFPT_PW_Data::set_phon_freq(int q_idx, const std::vector<double>& freq)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_phon_freq");
    ModuleBase::timer::start("DFPT_PW_Data", "set_phon_freq");
    if (q_idx >= static_cast<int>(phon_freq_.size()))
    {
        phon_freq_.resize(q_idx + 1);
    }
    phon_freq_[q_idx] = freq;
    ModuleBase::timer::end("DFPT_PW_Data", "set_phon_freq");
}

std::vector<double> DFPT_PW_Data::get_phon_freq(int q_idx) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_phon_freq");
    ModuleBase::timer::start("DFPT_PW_Data", "get_phon_freq");
    if (q_idx < static_cast<int>(phon_freq_.size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_phon_freq");
        return phon_freq_[q_idx];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_phon_freq");
    return std::vector<double>();
}

void DFPT_PW_Data::set_loto_dir(const ModuleBase::Vector3<double>& dir)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_loto_dir");
    ModuleBase::timer::start("DFPT_PW_Data", "set_loto_dir");
    const double norm = std::sqrt(dir * dir);
    const double null_norm_tol = 1.0e-10; ///< empirical parameter: norm below which the direction input is null
    if (norm < null_norm_tol)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_loto_dir");
        return; // keep the current direction on a null input
    }
    loto_dir_ = dir / norm;
    ModuleBase::timer::end("DFPT_PW_Data", "set_loto_dir");
}

void DFPT_PW_Data::set_dielectric(const ModuleBase::matrix& eps)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dielectric");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dielectric");
    dielectric_ = eps;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dielectric");
}

ModuleBase::matrix DFPT_PW_Data::get_dielectric() const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dielectric");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dielectric");
    ModuleBase::timer::end("DFPT_PW_Data", "get_dielectric");
    return dielectric_;
}

void DFPT_PW_Data::set_born(int atom_idx, const ModuleBase::matrix& z)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_born");
    ModuleBase::timer::start("DFPT_PW_Data", "set_born");
    if (atom_idx >= static_cast<int>(born_.size()))
    {
        born_.resize(atom_idx + 1);
    }
    born_[atom_idx] = z;
    ModuleBase::timer::end("DFPT_PW_Data", "set_born");
}

ModuleBase::matrix DFPT_PW_Data::get_born(int atom_idx) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_born");
    ModuleBase::timer::start("DFPT_PW_Data", "get_born");
    if (atom_idx < static_cast<int>(born_.size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_born");
        return born_[atom_idx];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_born");
    return ModuleBase::matrix();
}

void DFPT_PW_Data::allocate_memory()
{
    ModuleBase::TITLE("DFPT_PW_Data", "allocate_memory");
    ModuleBase::timer::start("DFPT_PW_Data", "allocate_memory");
    dynmat_.resize(get_nq());
    phon_freq_.resize(get_nq());
    born_.resize(nat_);
    ModuleBase::timer::end("DFPT_PW_Data", "allocate_memory");
}

void DFPT_PW_Data::deallocate_memory()
{
    ModuleBase::TITLE("DFPT_PW_Data", "deallocate_memory");
    ModuleBase::timer::start("DFPT_PW_Data", "deallocate_memory");
    dynmat_.clear();
    phon_freq_.clear();
    born_.clear();
    docc_.clear();
    dv_r_.clear();
    dv_recip_c_.clear();
    dv_rc_.clear();
    dpsi_.clear();
    converged_.clear();
    residuals_.clear();
    current_iter_.clear();
    ModuleBase::timer::end("DFPT_PW_Data", "deallocate_memory");
}

} // namespace ModuleDFPT
