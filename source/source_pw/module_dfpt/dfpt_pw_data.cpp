// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#include "dfpt_pw_data.h"
#include "source_lcao/module_dftu/dftu.h"

namespace ModuleDFPT {

DFPT_PW_Data::DFPT_PW_Data() {}

DFPT_PW_Data::~DFPT_PW_Data() {
    clean();
}

void DFPT_PW_Data::init(ModuleCell::QList* qlist, int nk, int nbands, int npw_max, 
                        int nrxx, int nspin, int nat, const Plus_U* dftu) {
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
}

void DFPT_PW_Data::clean() {
    deallocate_memory();
    is_initialized_ = false;
}

bool DFPT_PW_Data::u_active() const {
    // locale initialization requires the LCAO orbital files; a pure-PW run
    // without them has dftu != nullptr (wired upstream) but is not usable.
    return with_u() && dftu_->is_locale_initialized();
}

void DFPT_PW_Data::set_docc(int q_idx, const std::vector<std::complex<double>>& occ) {
    if (q_idx < 0) {
        return;
    }
    if (q_idx >= static_cast<int>(docc_.size())) {
        docc_.resize(q_idx + 1);
    }
    docc_[q_idx] = occ;
}

void DFPT_PW_Data::set_vsc_r(int atom_idx, int dir,
                             const std::vector<std::complex<double>>& v) {
    if (atom_idx < 0 || dir < 0 || dir >= 3) {
        return;
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot >= vsc_r_.size()) {
        vsc_r_.resize(slot + 1);
    }
    vsc_r_[slot] = v;
}

std::vector<std::complex<double>> DFPT_PW_Data::get_vsc_r(int atom_idx, int dir) const {
    if (atom_idx < 0 || dir < 0 || dir >= 3) {
        return std::vector<std::complex<double>>();
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot < vsc_r_.size()) {
        return vsc_r_[slot];
    }
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_dpsi_disp(
    int atom_idx, int dir,
    const std::vector<std::vector<std::vector<std::complex<double>>>>& d) {
    if (atom_idx < 0 || dir < 0 || dir >= 3) {
        return;
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot >= dpsi_disp_.size()) {
        dpsi_disp_.resize(slot + 1);
    }
    dpsi_disp_[slot] = d;
}

std::vector<std::vector<std::vector<std::complex<double>>>>
DFPT_PW_Data::get_dpsi_disp(int atom_idx, int dir) const {
    if (atom_idx < 0 || dir < 0 || dir >= 3) {
        return std::vector<std::vector<std::vector<std::complex<double>>>>();
    }
    const size_t slot = static_cast<size_t>(3 * atom_idx + dir);
    if (slot < dpsi_disp_.size()) {
        return dpsi_disp_[slot];
    }
    return std::vector<std::vector<std::vector<std::complex<double>>>>();
}

std::vector<std::complex<double>> DFPT_PW_Data::get_docc(int q_idx) const {
    if (q_idx >= 0 && q_idx < static_cast<int>(docc_.size())) {
        return docc_[q_idx];
    }
    return std::vector<std::complex<double>>();
}

int DFPT_PW_Data::get_nq() const {
    return qlist_->get_nq();
}

ModuleBase::Vector3<double> DFPT_PW_Data::get_qvec(int q_idx) const {
    return qlist_->get_q(q_idx);
}

int DFPT_PW_Data::get_nirr(int q_idx) const {
    return qlist_->get_nirr(q_idx);
}

std::vector<int> DFPT_PW_Data::get_irrep_modes(int q_idx, int irrep) const {
    return qlist_->get_irrep_modes(q_idx, irrep);
}

void DFPT_PW_Data::set_dpsi(int q_idx, int k_idx, int band_idx, 
                             const std::vector<std::complex<double>>& psi) {
    if (q_idx < 0 || k_idx < 0 || band_idx < 0) {
        return;
    }
    if (q_idx >= static_cast<int>(dpsi_.size())) {
        dpsi_.resize(q_idx + 1);
    }
    if (k_idx >= static_cast<int>(dpsi_[q_idx].size())) {
        dpsi_[q_idx].resize(k_idx + 1);
    }
    if (band_idx >= static_cast<int>(dpsi_[q_idx][k_idx].size())) {
        dpsi_[q_idx][k_idx].resize(band_idx + 1);
    }
    dpsi_[q_idx][k_idx][band_idx] = psi;
}

std::vector<std::complex<double>> DFPT_PW_Data::get_dpsi(int q_idx, int k_idx, int band_idx) const {
    if (q_idx >= 0 && k_idx >= 0 && band_idx >= 0 &&
        q_idx < static_cast<int>(dpsi_.size()) &&
        k_idx < static_cast<int>(dpsi_[q_idx].size()) &&
        band_idx < static_cast<int>(dpsi_[q_idx][k_idx].size()))
    {
        return dpsi_[q_idx][k_idx][band_idx];
    }
    return std::vector<std::complex<double>>();
}

psi::Psi<std::complex<double>>& DFPT_PW_Data::get_dpsi_obj(int q_idx) {
    static psi::Psi<std::complex<double>> dummy;
    (void)q_idx;
    return dummy;
}

void DFPT_PW_Data::set_drho_r(int q_idx, int spin, const std::vector<double>& rho) {
    if (q_idx < 0 || spin < 0) { return; }
    if (q_idx >= static_cast<int>(drho_r_.size())) {
        drho_r_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(drho_r_[q_idx].size())) {
        drho_r_[q_idx].resize(spin + 1);
    }
    drho_r_[q_idx][spin] = rho;
}

std::vector<double> DFPT_PW_Data::get_drho_r(int q_idx, int spin) const {
    if (q_idx >= 0 && spin >= 0 &&
        q_idx < static_cast<int>(drho_r_.size()) &&
        spin < static_cast<int>(drho_r_[q_idx].size()))
    {
        return drho_r_[q_idx][spin];
    }
    return std::vector<double>();
}

void DFPT_PW_Data::set_drho_g(int q_idx, int spin, const std::vector<std::complex<double>>& rho) {
    if (q_idx < 0 || spin < 0) { return; }
    if (q_idx >= static_cast<int>(drho_g_.size())) {
        drho_g_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(drho_g_[q_idx].size())) {
        drho_g_[q_idx].resize(spin + 1);
    }
    drho_g_[q_idx][spin] = rho;
}

std::vector<std::complex<double>> DFPT_PW_Data::get_drho_g(int q_idx, int spin) const {
    if (q_idx >= 0 && spin >= 0 &&
        q_idx < static_cast<int>(drho_g_.size()) &&
        spin < static_cast<int>(drho_g_[q_idx].size()))
    {
        return drho_g_[q_idx][spin];
    }
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_dv_r(int q_idx, int spin, const std::vector<double>& v) {
    if (q_idx < 0 || spin < 0) { return; }
    if (q_idx >= static_cast<int>(dv_r_.size())) {
        dv_r_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(dv_r_[q_idx].size())) {
        dv_r_[q_idx].resize(spin + 1);
    }
    dv_r_[q_idx][spin] = v;
}

std::vector<double> DFPT_PW_Data::get_dv_r(int q_idx, int spin) const {
    if (q_idx >= 0 && spin >= 0 &&
        q_idx < static_cast<int>(dv_r_.size()) &&
        spin < static_cast<int>(dv_r_[q_idx].size()))
    {
        return dv_r_[q_idx][spin];
    }
    return std::vector<double>();
}

void DFPT_PW_Data::set_dv_recip_c(int q_idx, int spin, const std::vector<std::complex<double>>& v) {
    if (q_idx < 0 || spin < 0) { return; }
    if (q_idx >= static_cast<int>(dv_recip_c_.size())) {
        dv_recip_c_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(dv_recip_c_[q_idx].size())) {
        dv_recip_c_[q_idx].resize(spin + 1);
    }
    dv_recip_c_[q_idx][spin] = v;
}

std::vector<std::complex<double>> DFPT_PW_Data::get_dv_recip_c(int q_idx, int spin) const {
    if (q_idx >= 0 && spin >= 0 &&
        q_idx < static_cast<int>(dv_recip_c_.size()) &&
        spin < static_cast<int>(dv_recip_c_[q_idx].size()))
    {
        return dv_recip_c_[q_idx][spin];
    }
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_dv_rc(int q_idx, int spin, const std::vector<std::complex<double>>& v) {
    if (q_idx < 0 || spin < 0) { return; }
    if (q_idx >= static_cast<int>(dv_rc_.size())) {
        dv_rc_.resize(q_idx + 1);
    }
    if (spin >= static_cast<int>(dv_rc_[q_idx].size())) {
        dv_rc_[q_idx].resize(spin + 1);
    }
    dv_rc_[q_idx][spin] = v;
}

std::vector<std::complex<double>> DFPT_PW_Data::get_dv_rc(int q_idx, int spin) const {
    if (q_idx >= 0 && spin >= 0 &&
        q_idx < static_cast<int>(dv_rc_.size()) &&
        spin < static_cast<int>(dv_rc_[q_idx].size()))
    {
        return dv_rc_[q_idx][spin];
    }
    return std::vector<std::complex<double>>();
}

void DFPT_PW_Data::set_dynmat(int q_idx, const ModuleBase::ComplexMatrix& dm) {
    if (q_idx >= static_cast<int>(dynmat_.size())) {
        dynmat_.resize(q_idx + 1);
    }
    dynmat_[q_idx] = dm;
}

ModuleBase::ComplexMatrix DFPT_PW_Data::get_dynmat(int q_idx) const {
    if (q_idx < static_cast<int>(dynmat_.size())) {
        return dynmat_[q_idx];
    }
    return ModuleBase::ComplexMatrix();
}

void DFPT_PW_Data::set_phon_freq(int q_idx, const std::vector<double>& freq) {
    if (q_idx >= static_cast<int>(phon_freq_.size())) {
        phon_freq_.resize(q_idx + 1);
    }
    phon_freq_[q_idx] = freq;
}

std::vector<double> DFPT_PW_Data::get_phon_freq(int q_idx) const {
    if (q_idx < static_cast<int>(phon_freq_.size())) {
        return phon_freq_[q_idx];
    }
    return std::vector<double>();
}

void DFPT_PW_Data::set_dielectric(const ModuleBase::matrix& eps) {
    dielectric_ = eps;
}

ModuleBase::matrix DFPT_PW_Data::get_dielectric() const {
    return dielectric_;
}

void DFPT_PW_Data::set_born(int atom_idx, const ModuleBase::matrix& z) {
    if (atom_idx >= static_cast<int>(born_.size())) {
        born_.resize(atom_idx + 1);
    }
    born_[atom_idx] = z;
}

ModuleBase::matrix DFPT_PW_Data::get_born(int atom_idx) const {
    if (atom_idx < static_cast<int>(born_.size())) {
        return born_[atom_idx];
    }
    return ModuleBase::matrix();
}

void DFPT_PW_Data::allocate_memory() {
    dynmat_.resize(get_nq());
    phon_freq_.resize(get_nq());
    born_.resize(nat_);
}

void DFPT_PW_Data::deallocate_memory() {
    dynmat_.clear();
    phon_freq_.clear();
    born_.clear();
    docc_.clear();
    dv_r_.clear();
    dv_recip_c_.clear();
    dv_rc_.clear();
    dpsi_.clear();
    residuals_.clear();
}

} // namespace ModuleDFPT