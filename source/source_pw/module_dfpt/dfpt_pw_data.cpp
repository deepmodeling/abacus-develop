#include "dfpt_pw_data.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include <cmath>
#include <sstream>

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

int DFPT_PW_Data::get_nq() const
{
    return qlist_->get_nq();
}

void validate_q_ordering(const ModuleCell::QList& qlist, bool compute_q0)
{
    ModuleBase::TITLE("ModuleDFPT", "validate_q_ordering");
    ModuleBase::timer::start("ModuleDFPT", "validate_q_ordering");
    if (!compute_q0)
    {
        ModuleBase::timer::end("ModuleDFPT", "validate_q_ordering");
        return;
    }
    const double gamma_tol = 1.0e-8; ///< empirical parameter: fractional-q Gamma tolerance (dfpt_phon.cpp value)
    const int nq = qlist.get_nq();
    const ModuleBase::Vector3<double> q_first = qlist.get_q(0);
    const bool first_is_gamma = std::abs(q_first.x) <= gamma_tol && std::abs(q_first.y) <= gamma_tol
                                && std::abs(q_first.z) <= gamma_tol;
    if (!first_is_gamma)
    {
        ModuleBase::timer::end("ModuleDFPT", "validate_q_ordering");
        ModuleBase::WARNING_QUIT("validate_q_ordering",
                                 "dfpt_compute_q0 requires Gamma (0 0 0) as the FIRST q point: "
                                 "the E-field response, dielectric tensor and Born charges are "
                                 "computed at q index 0; put Gamma first in the q-point file.");
    }
    for (int iq = 1; iq < nq; ++iq)
    {
        const ModuleBase::Vector3<double> qi = qlist.get_q(iq);
        const bool is_gamma = std::abs(qi.x) <= gamma_tol && std::abs(qi.y) <= gamma_tol
                              && std::abs(qi.z) <= gamma_tol;
        if (is_gamma)
        {
            std::stringstream msg;
            msg << "duplicate Gamma q point at index " << iq
                << " in the q-point list: with dfpt_compute_q0 the polar tensors and LO-TO "
                   "correction attach to the first Gamma entry only; remove the duplicate "
                   "Gamma point.";
            ModuleBase::timer::end("ModuleDFPT", "validate_q_ordering");
            ModuleBase::WARNING_QUIT("validate_q_ordering", msg.str());
        }
    }
    ModuleBase::timer::end("ModuleDFPT", "validate_q_ordering");
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

void DFPT_PW_Data::set_dynmat_raw(int q_idx, const ModuleBase::ComplexMatrix& dm)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_dynmat_raw");
    ModuleBase::timer::start("DFPT_PW_Data", "set_dynmat_raw");
    if (q_idx >= static_cast<int>(dynmat_raw_.size()))
    {
        const int raw_size = q_idx + 1;
        dynmat_raw_.resize(raw_size);
    }
    dynmat_raw_[q_idx] = dm;
    ModuleBase::timer::end("DFPT_PW_Data", "set_dynmat_raw");
}

ModuleBase::ComplexMatrix DFPT_PW_Data::get_dynmat_raw(int q_idx) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_dynmat_raw");
    ModuleBase::timer::start("DFPT_PW_Data", "get_dynmat_raw");
    if (q_idx < static_cast<int>(dynmat_raw_.size()) && dynmat_raw_[q_idx].nr > 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_dynmat_raw");
        return dynmat_raw_[q_idx];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_dynmat_raw");
    return get_dynmat(q_idx); // fallback: the working copy
}

bool DFPT_PW_Data::has_dynmat_raw(int q_idx) const
{
    return q_idx >= 0 && q_idx < static_cast<int>(dynmat_raw_.size()) && dynmat_raw_[q_idx].nr > 0;
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
    dynmat_raw_.clear();
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

void DFPT_PW_Data::set_drhoc_g(int q_idx, const std::vector<std::complex<double>>& rho_c)
{
    ModuleBase::TITLE("DFPT_PW_Data", "set_drhoc_g");
    ModuleBase::timer::start("DFPT_PW_Data", "set_drhoc_g");
    if (q_idx < 0)
    {
        ModuleBase::timer::end("DFPT_PW_Data", "set_drhoc_g");
        return;
    }
    if (q_idx >= static_cast<int>(drhoc_g_.size()))
    {
        drhoc_g_.resize(q_idx + 1);
    }
    drhoc_g_[q_idx] = rho_c;
    ModuleBase::timer::end("DFPT_PW_Data", "set_drhoc_g");
}

std::vector<std::complex<double>> DFPT_PW_Data::get_drhoc_g(int q_idx) const
{
    ModuleBase::TITLE("DFPT_PW_Data", "get_drhoc_g");
    ModuleBase::timer::start("DFPT_PW_Data", "get_drhoc_g");
    if (q_idx >= 0 && q_idx < static_cast<int>(drhoc_g_.size()))
    {
        ModuleBase::timer::end("DFPT_PW_Data", "get_drhoc_g");
        return drhoc_g_[q_idx];
    }
    ModuleBase::timer::end("DFPT_PW_Data", "get_drhoc_g");
    return {};
}

} // namespace ModuleDFPT
