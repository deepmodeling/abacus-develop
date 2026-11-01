// ============================================================
// Response-field slots of DFPT_PW_Data: the bound-safe getters/setters
// of the per-direction / per-q response fields -- u_active, docc,
// v_sc, the displacement / E-field response stashes (drho_disp,
// dpsi_disp, pos_resp, dpsi_efield) and the drho / dv slot family.
// The driver state (lifecycle, q delegation, dpsi ledger, dynmat /
// phonon tensors) stays in dfpt_pw_data.cpp.
// ============================================================

#include "dfpt_pw_data.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <complex>
#include <vector>

namespace ModuleDFPT
{


bool DFPT_PW_Data::u_active() const
{
    ModuleBase::TITLE("DFPT_PW_Data", "u_active");
    ModuleBase::timer::start("DFPT_PW_Data", "u_active");
    // a usable provider has its occupation matrices initialized (the ground
    // state does this when DFT+U actually runs); a wired provider without
    // them (e.g. a default-constructed reservation) stays inactive.
    ModuleBase::timer::end("DFPT_PW_Data", "u_active");
    return with_u() && dftu_->is_occmat_ready();
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

void DFPT_PW_Data::set_drho_disp(int atom_idx, int dir, const std::vector<std::complex<double>>& rho)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_drho_disp");
    ModuleBase::timer::start("DFPT_PW_Data", "set_drho_disp");
    if (atom_idx < 0 || dir < 0 || dir >= 3)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_drho_disp");
        return;
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot >= drho_disp_.size())
    {
        drho_disp_.resize(slot + 1);
    }
    drho_disp_[slot] = rho;
    ModuleBase::timer::end("DFPT_PW_Data", "set_drho_disp");
}

std::vector<std::complex<double>> DFPT_PW_Data::get_drho_disp(int atom_idx, int dir) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_drho_disp");
    ModuleBase::timer::start("DFPT_PW_Data", "get_drho_disp");
    if (atom_idx < 0 || dir < 0 || dir >= 3)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_drho_disp");
        return std::vector<std::complex<double>>();
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot < drho_disp_.size())
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_drho_disp");
        return drho_disp_[slot];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_drho_disp");
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


} // namespace ModuleDFPT
