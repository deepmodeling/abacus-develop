/**
 * @file diagonalization_engine.cpp
 * @brief Implementations of DiagonalizationEngine strategies for DeltaSpin.
 *
 * @par FullSpaceDiagonalizer
 * Direct call to HSolverLCAO::solve(). Ground truth, always correct.
 *
 * @par SubspaceDiagonalizer
 * Builds H_sub = H0_sub + dl * P_I_sub, diagonalizes in subspace,
 * rotates wavefunctions, computes Mi from full-space DMR.
 *
 * @par FirstOrderResponseEngine
 * Shifts eigenvalues: eps_new = eps_old + dl * diag(P_I_sub).
 * Wavefunctions frozen. Fastest but least accurate.
 */

#include "diagonalization_engine.h"
#include "spin_constrain.h"

#include <algorithm>
#include <cstring>
#include <iostream>

#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#ifdef __LCAO
#include "source_hsolver/hsolver_lcao.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/kernels/hegvd_op.h"
#include "source_estate/elecstate_tools.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#endif

namespace spinconstrain
{

// ===================================================================
// SubspaceCache implementation
// ===================================================================

SubspaceCache::SubspaceCache() = default;

void SubspaceCache::build(
    int nk,
    int nbands,
    int nat,
    const std::complex<double>* H0_sub_raw,
    const std::complex<double>* S_sub_raw,
    std::vector<std::vector<std::vector<std::complex<double>>>> P_I_sub_all,
    const std::vector<double>& ekb_ref_all,
    const std::vector<ModuleBase::Vector3<double>>& lambda_ref)
{
    nk_ = nk;
    nbands_ = nbands;
    nat_ = nat;

    const int nn = nbands * nbands;

    // Copy H0_sub and S_sub
    H0_sub_.assign(H0_sub_raw, H0_sub_raw + nk * nn);
    S_sub_.assign(S_sub_raw, S_sub_raw + nk * nn);

    // P_I_sub: deep copy
    P_I_sub_ = std::move(P_I_sub_all);

    // Eigenvalues
    ekb_ref_ = ekb_ref_all;

    // Lambda reference
    lambda_ref_ = lambda_ref;

    valid_ = true;
}

void SubspaceCache::clear()
{
    nk_ = 0;
    nbands_ = 0;
    nat_ = 0;
    valid_ = false;

    H0_sub_.clear();
    S_sub_.clear();
    P_I_sub_.clear();
    ekb_ref_.clear();
    lambda_ref_.clear();
}

const std::complex<double>* SubspaceCache::H0_sub(int ik) const
{
    if (!valid_ || ik < 0 || ik >= nk_) return nullptr;
    return H0_sub_.data() + ik * nbands_ * nbands_;
}

const std::complex<double>* SubspaceCache::S_sub(int ik) const
{
    if (!valid_ || ik < 0 || ik >= nk_) return nullptr;
    return S_sub_.data() + ik * nbands_ * nbands_;
}

const std::complex<double>* SubspaceCache::P_I_sub(int ik, int iat) const
{
    if (!valid_ || ik < 0 || ik >= nk_ || iat < 0 || iat >= nat_) return nullptr;
    if (P_I_sub_[ik][iat].empty()) return nullptr;
    return P_I_sub_[ik][iat].data();
}

const double* SubspaceCache::ekb_ref(int ik) const
{
    if (!valid_ || ik < 0 || ik >= nk_) return nullptr;
    return ekb_ref_.data() + ik * nbands_;
}

const std::vector<ModuleBase::Vector3<double>>& SubspaceCache::lambda_ref() const
{
    return lambda_ref_;
}

// ===================================================================
// FullSpaceDiagonalizer implementation
// ===================================================================

FullSpaceDiagonalizer::FullSpaceDiagonalizer(SpinConstrain<std::complex<double>>& sc)
    : sc_(sc)
{
}

DiagonalizationResult FullSpaceDiagonalizer::solve(int i_step)
{
    ModuleBase::timer::start("FullSpaceDiagonalizer", "solve");

    DiagonalizationResult result{};

#ifdef __LCAO
    if (sc_.get_basis_type() != "lcao")
    {
        // FullSpaceDiagonalizer is LCAO-only; signal failure
        result.success = false;
        ModuleBase::timer::end("FullSpaceDiagonalizer", "solve");
        return result;
    }

    psi::Psi<std::complex<double>>* psi_t
        = static_cast<psi::Psi<std::complex<double>>*>(sc_.psi);
    hamilt::Hamilt<std::complex<double>>* hamilt_t
        = static_cast<hamilt::Hamilt<std::complex<double>>*>(sc_.p_hamilt);

    // Update DeltaSpin operator with current lambda
    if (sc_.get_nspin() == 2)
    {
        dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(
            sc_.p_operator)
            ->update_lambda();
    }
    else if (sc_.get_nspin() == 4)
    {
        dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>* >(
            sc_.p_operator)
            ->update_lambda();
    }

    // Full diagonalization without charge update (last param = true)
    hsolver::HSolverLCAO<std::complex<double>> hsolver_t(sc_.ParaV, sc_.get_ks_solver());
    hsolver_t.solve(hamilt_t, *psi_t, sc_.pelec, *sc_.dm_, *sc_.pelec->charge,
                    sc_.get_nspin(), true);

    // Compute weights from new eigenvalues
    elecstate::calculate_weights(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->klist,
                                 sc_.pelec->eferm, sc_.pelec->f_en,
                                 sc_.pelec->nelec_spin, sc_.pelec->skip_weights);
    elecstate::calEBand(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->f_en);

    // Compute magnetic moments
    sc_.cal_mi_lcao(i_step);

    result.success = true;
    result.used_subspace_approximation = false;
    result.nbands = sc_.get_nbands();
    result.nk = psi_t->get_nk();
#else
    result.success = false;
#endif

    ModuleBase::timer::end("FullSpaceDiagonalizer", "solve");
    return result;
}

bool FullSpaceDiagonalizer::build_subspace(
    const std::vector<ModuleBase::Vector3<double>>& lambda_ref)
{
    // FullSpaceDiagonalizer does not use subspace cache
    (void)lambda_ref;
    return false;
}

bool FullSpaceDiagonalizer::has_subspace() const
{
    return false;
}

void FullSpaceDiagonalizer::clear_subspace()
{
    // No-op: FullSpaceDiagonalizer has no subspace cache
}

// ===================================================================
// SubspaceDiagonalizer implementation
// ===================================================================

SubspaceDiagonalizer::SubspaceDiagonalizer(SpinConstrain<std::complex<double>>& sc)
    : sc_(sc)
{
}

DiagonalizationResult SubspaceDiagonalizer::solve(int i_step)
{
    ModuleBase::timer::start("SubspaceDiagonalizer", "solve");

    DiagonalizationResult result{};

#ifdef __LCAO
    if (sc_.get_basis_type() != "lcao")
    {
        result.success = false;
        ModuleBase::timer::end("SubspaceDiagonalizer", "solve");
        return result;
    }

    if (!cache_.is_valid())
    {
        result.success = false;
        ModuleBase::timer::end("SubspaceDiagonalizer", "solve");
        return result;
    }

    const int nk = cache_.nk();
    const int nbands = cache_.nbands();
    const int nn = nbands * nbands;
    const int nat = sc_.get_nat();

    // Resize temporary buffers if needed
    if (static_cast<int>(h_tmp_.size()) < nn)
    {
        h_tmp_.resize(nn);
        s_tmp_.resize(nn);
        s_copy_.resize(nn);
        eigenvalues_.resize(nbands);
        eigenvectors_.resize(nn);
    }

    // Update DeltaSpin operator with current lambda
    if (sc_.get_nspin() == 2)
    {
        dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(
            sc_.p_operator)
            ->update_lambda();
    }
    else if (sc_.get_nspin() == 4)
    {
        dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>* >(
            sc_.p_operator)
            ->update_lambda();
    }

    // Storage for eigenvectors per k-point
    std::vector<std::vector<std::complex<double>>> vcc_all(nk);

    double max_eig_change = 0.0;

    for (int ik = 0; ik < nk; ik++)
    {
        // Copy H0_sub from cache to temporary
        std::memcpy(h_tmp_.data(), cache_.H0_sub(ik), sizeof(std::complex<double>) * nn);

        // Build PI_sub for this k-point from cache
        std::vector<std::vector<std::complex<double>>> PI_sub_for_ik(nat);
        for (int iat = 0; iat < nat; iat++)
        {
            const std::complex<double>* p = cache_.P_I_sub(ik, iat);
            if (p != nullptr)
            {
                PI_sub_for_ik[iat].assign(p, p + nn);
            }
        }

        sc_.calculate_delta_hcc_lcao(h_tmp_.data(), PI_sub_for_ik,
                                     sc_.get_sc_lambda().data(), nbands, ik, true);

        // Copy S_sub for diag_hegvd (it modifies S in place)
        std::memcpy(s_copy_.data(), cache_.S_sub(ik), sizeof(std::complex<double>) * nn);

        // Diagonalize: H_sub * V = S_sub * V * eps
        hsolver::DiagoIterAssist<std::complex<double>>::diag_hegvd(
            nbands, nbands, h_tmp_.data(), s_copy_.data(), nbands,
            eigenvalues_.data(), eigenvectors_.data());

        // Store eigenvectors
        vcc_all[ik].assign(eigenvectors_.data(), eigenvectors_.data() + nn);

        // Update eigenvalues in pelec
        for (int ib = 0; ib < nbands; ib++)
        {
            double old_eig = cache_.ekb_ref(ik)[ib];
            sc_.pelec->ekb(ik, ib) = eigenvalues_[ib];
            double change = std::abs(eigenvalues_[ib] - old_eig);
            if (change > max_eig_change) max_eig_change = change;
        }
    }

    // Compute weights from new eigenvalues
    elecstate::calculate_weights(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->klist,
                                 sc_.pelec->eferm, sc_.pelec->f_en,
                                 sc_.pelec->nelec_spin, sc_.pelec->skip_weights);
    elecstate::calEBand(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->f_en);

    // Compute Mi by rotating psi, computing full-space DMR, then restoring psi.
    sc_.cal_mi_lcao_subspace(vcc_all, nbands, nk, sc_.get_npol());

    result.success = true;
    result.used_subspace_approximation = true;
    result.max_eigenvalue_change = max_eig_change;
    result.nbands = nbands;
    result.nk = nk;
#else
    result.success = false;
#endif

    ModuleBase::timer::end("SubspaceDiagonalizer", "solve");
    return result;
}

bool SubspaceDiagonalizer::build_subspace(
    const std::vector<ModuleBase::Vector3<double>>& lambda_ref)
{
#ifdef __LCAO
    if (sc_.get_basis_type() != "lcao") return false;

    const int nk = sc_.pelec->klist->get_nks();
    const int nbands = sc_.get_nbands();
    const int nat = sc_.get_nat();
    const int nlocal = sc_.ParaV->nrow;
    const int nn = nbands * nbands;

    psi::Psi<std::complex<double>>* psi_t
        = static_cast<psi::Psi<std::complex<double>>*>(sc_.psi);
    hamilt::Hamilt<std::complex<double>>* hamilt_t
        = static_cast<hamilt::Hamilt<std::complex<double>>*>(sc_.p_hamilt);

    // CRITICAL: Run full HSolverLCAO diagonalization BEFORE building subspace.
    // This updates psi to be the eigenvectors of H(lambda_ref). The subspace
    // C^† H C is only diagonal (and the perturbation approximation valid)
    // when C contains the eigenvectors of H(lambda_ref).
    // Without this step, psi is stale (from BFGS trial step) and H0_sub
    // will have large off-diagonal elements, breaking the perturbation theory.
    {
        if (sc_.get_nspin() == 2)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double> > *>(
                sc_.p_operator)->update_lambda();
        }
        else if (sc_.get_nspin() == 4)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double> > > *>(
                sc_.p_operator)->update_lambda();
        }

        hsolver::HSolverLCAO<std::complex<double>> hsolver_t(sc_.ParaV, sc_.get_ks_solver());
        hsolver_t.solve(hamilt_t, *psi_t, sc_.pelec, *sc_.dm_, *sc_.pelec->charge,
                        sc_.get_nspin(), true);

        elecstate::calculate_weights(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->klist,
                                     sc_.pelec->eferm, sc_.pelec->f_en,
                                     sc_.pelec->nelec_spin, sc_.pelec->skip_weights);
        elecstate::calEBand(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->f_en);
    }

    // Clear existing cache
    cache_.clear();

    // Allocate raw arrays for building
    std::vector<std::complex<double>> H0_sub_raw(nk * nn);
    std::vector<std::complex<double>> S_sub_raw(nk * nn);
    std::vector<std::vector<std::vector<std::complex<double>>>> P_I_sub_all(nk);

    // Compute H0_sub, S_sub, P_I_sub for each k-point
    for (int ik = 0; ik < nk; ik++)
    {
        psi_t->fix_k(ik);

        // H0_sub = C^dag H C, S_sub = C^dag S C
        // NOW C contains eigenvectors of H(lambda_ref), so H0_sub is diagonal
        sc_.calculate_lcao_sub_hs(
            sc_.p_hamilt, *psi_t, sc_.ParaV,
            H0_sub_raw.data() + ik * nn,
            S_sub_raw.data() + ik * nn,
            ik, nbands, nlocal);

        // P_I_sub = C^dag D_I C for each constrained atom
        if (sc_.get_nspin() == 2)
        {
            auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double> > *>(
                sc_.p_operator);
            dspin_op->cal_PI_sub(
                sc_.kv_.kvec_d[ik],
                psi_t->get_pointer(),
                nbands,
                P_I_sub_all[ik]);
        }
        else if (sc_.get_nspin() == 4)
        {
            auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>> > *>(
                sc_.p_operator);
            dspin_op->cal_PI_sub(
                sc_.kv_.kvec_d[ik],
                psi_t->get_pointer(),
                nbands,
                P_I_sub_all[ik]);
        }
    }

    // Collect eigenvalues
    std::vector<double> ekb_ref_all(nk * nbands);
    for (int ik = 0; ik < nk; ik++)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            ekb_ref_all[ik * nbands + ib] = sc_.pelec->ekb(ik, ib);
        }
    }

    // Build the cache
    cache_.build(nk, nbands, nat,
                 H0_sub_raw.data(), S_sub_raw.data(),
                 std::move(P_I_sub_all),
                 ekb_ref_all, lambda_ref);

    // CRITICAL: Synchronize lcao_lambda_in_sub_ so that calculate_delta_hcc_lcao
    // (which uses this->lcao_lambda_in_sub_ internally) computes the correct
    // delta = lambda - lambda_ref. Without this, calculate_delta_hcc_lcao would
    // use a stale or uninitialized lcao_lambda_in_sub_, producing a completely
    // wrong H_sub and catastrophic Mi errors (RMS jumps from ~0.003 to ~0.35).
    sc_.lcao_lambda_in_sub_ = lambda_ref;

    // Compute Mi from the full diagonalization result (consistent with original i_step==-2)
    sc_.cal_mi_lcao(-2);

    return true;
#else
    (void)lambda_ref;
    return false;
#endif
}

bool SubspaceDiagonalizer::has_subspace() const
{
    return cache_.is_valid();
}

void SubspaceDiagonalizer::clear_subspace()
{
    cache_.clear();
}

// ===================================================================
// FirstOrderResponseEngine implementation
// ===================================================================

FirstOrderResponseEngine::FirstOrderResponseEngine(SpinConstrain<std::complex<double>>& sc)
    : sc_(sc)
{
}

DiagonalizationResult FirstOrderResponseEngine::solve(int i_step)
{
    ModuleBase::timer::start("FirstOrderResponseEngine", "solve");

    DiagonalizationResult result{};

#ifdef __LCAO
    if (sc_.get_basis_type() != "lcao")
    {
        result.success = false;
        ModuleBase::timer::end("FirstOrderResponseEngine", "solve");
        return result;
    }

    if (!cache_.is_valid())
    {
        result.success = false;
        ModuleBase::timer::end("FirstOrderResponseEngine", "solve");
        return result;
    }

    const int nk = cache_.nk();
    const int nbands = cache_.nbands();
    const int nat = sc_.get_nat();
    const int nn = nbands * nbands;

    const auto& lambda = sc_.get_sc_lambda();
    const auto& lambda_ref = cache_.lambda_ref();

    double max_eig_change = 0.0;

    // First-order eigenvalue shift: delta_eps = sum_I dl_I * diag(P_I_sub)_ib
    for (int ik = 0; ik < nk; ik++)
    {
        int spin_sign = sc_.get_spin_sign(ik);

        for (int ib = 0; ib < nbands; ib++)
        {
            double delta_epsilon = 0.0;

            for (int iat = 0; iat < nat; iat++)
            {
                const std::complex<double>* p = cache_.P_I_sub(ik, iat);
                if (p == nullptr) continue;

                double p_diag = p[ib + ib * nbands].real();
                double dl = lambda[iat][2] - lambda_ref[iat][2];
                delta_epsilon += dl * p_diag;
            }

            double old_eig = cache_.ekb_ref(ik)[ib];
            double new_eig = old_eig - spin_sign * delta_epsilon;
            sc_.pelec->ekb(ik, ib) = new_eig;

            double change = std::abs(new_eig - old_eig);
            if (change > max_eig_change) max_eig_change = change;
        }
    }

    // Compute weights from shifted eigenvalues
    elecstate::calculate_weights(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->klist,
                                 sc_.pelec->eferm, sc_.pelec->f_en,
                                 sc_.pelec->nelec_spin, sc_.pelec->skip_weights);
    elecstate::calEBand(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->f_en);

    // Compute Mi from subspace: first-order uses V = I (no rotation), so
    // psi unchanged, DMK/DMR computed from current weights.
    // Use cal_mi_lcao_subspace with identity eigenvectors.
    std::vector<std::vector<std::complex<double>>> vcc_identity(nk,
        std::vector<std::complex<double>>(nn, {0.0, 0.0}));
    for (int ik = 0; ik < nk; ik++)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            vcc_identity[ik][ib + ib * nbands] = {1.0, 0.0};
        }
    }

    sc_.cal_mi_lcao_subspace(vcc_identity, nbands, nk, sc_.get_npol());

    result.success = true;
    result.used_subspace_approximation = true;
    result.max_eigenvalue_change = max_eig_change;
    result.nbands = nbands;
    result.nk = nk;
#else
    result.success = false;
#endif

    ModuleBase::timer::end("FirstOrderResponseEngine", "solve");
    return result;
}

bool FirstOrderResponseEngine::build_subspace(
    const std::vector<ModuleBase::Vector3<double>>& lambda_ref)
{
#ifdef __LCAO
    if (sc_.get_basis_type() != "lcao") return false;

    const int nk = sc_.pelec->klist->get_nks();
    const int nbands = sc_.get_nbands();
    const int nat = sc_.get_nat();
    const int nlocal = sc_.ParaV->nrow;
    const int nn = nbands * nbands;

    psi::Psi<std::complex<double>>* psi_t
        = static_cast<psi::Psi<std::complex<double>>*>(sc_.psi);

    // CRITICAL: Same as SubspaceDiagonalizer - must run full diagonalization
    // first to update psi to eigenvectors of H(lambda_ref).
    {
        if (sc_.get_nspin() == 2)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double> > *>(
                sc_.p_operator)->update_lambda();
        }
        else if (sc_.get_nspin() == 4)
        {
            dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double> > > *>(
                sc_.p_operator)->update_lambda();
        }

        hamilt::Hamilt<std::complex<double>>* hamilt_t
            = static_cast<hamilt::Hamilt<std::complex<double>>*>(sc_.p_hamilt);
        hsolver::HSolverLCAO<std::complex<double>> hsolver_t(sc_.ParaV, sc_.get_ks_solver());
        hsolver_t.solve(hamilt_t, *psi_t, sc_.pelec, *sc_.dm_, *sc_.pelec->charge,
                        sc_.get_nspin(), true);

        elecstate::calculate_weights(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->klist,
                                     sc_.pelec->eferm, sc_.pelec->f_en,
                                     sc_.pelec->nelec_spin, sc_.pelec->skip_weights);
        elecstate::calEBand(sc_.pelec->ekb, sc_.pelec->wg, sc_.pelec->f_en);
    }

    cache_.clear();

    std::vector<std::complex<double>> H0_sub_raw(nk * nn);
    std::vector<std::complex<double>> S_sub_raw(nk * nn);
    std::vector<std::vector<std::vector<std::complex<double>>>> P_I_sub_all(nk);

    for (int ik = 0; ik < nk; ik++)
    {
        psi_t->fix_k(ik);

        sc_.calculate_lcao_sub_hs(
            sc_.p_hamilt, *psi_t, sc_.ParaV,
            H0_sub_raw.data() + ik * nn,
            S_sub_raw.data() + ik * nn,
            ik, nbands, nlocal);

        if (sc_.get_nspin() == 2)
        {
            auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double> > *>(
                sc_.p_operator);
            dspin_op->cal_PI_sub(
                sc_.kv_.kvec_d[ik],
                psi_t->get_pointer(),
                nbands,
                P_I_sub_all[ik]);
        }
        else if (sc_.get_nspin() == 4)
        {
            auto* dspin_op = dynamic_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>> > *>(
                sc_.p_operator);
            dspin_op->cal_PI_sub(
                sc_.kv_.kvec_d[ik],
                psi_t->get_pointer(),
                nbands,
                P_I_sub_all[ik]);
        }
    }

    std::vector<double> ekb_ref_all(nk * nbands);
    for (int ik = 0; ik < nk; ik++)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            ekb_ref_all[ik * nbands + ib] = sc_.pelec->ekb(ik, ib);
        }
    }

    cache_.build(nk, nbands, nat,
                 H0_sub_raw.data(), S_sub_raw.data(),
                 std::move(P_I_sub_all),
                 ekb_ref_all, lambda_ref);

    sc_.lcao_lambda_in_sub_ = lambda_ref;

    sc_.cal_mi_lcao(-2);

    return true;
#else
    (void)lambda_ref;
    return false;
#endif
}

bool FirstOrderResponseEngine::has_subspace() const
{
    return cache_.is_valid();
}

void FirstOrderResponseEngine::clear_subspace()
{
    cache_.clear();
}

// ===================================================================
// Factory function
// ===================================================================

std::unique_ptr<DiagonalizationEngine> create_diagonalization_engine(
    DiagonalizationStrategy strategy,
    SpinConstrain<std::complex<double>>& sc)
{
    switch (strategy)
    {
        case DiagonalizationStrategy::FullSpace:
        {
            std::unique_ptr<DiagonalizationEngine> ptr(new FullSpaceDiagonalizer(sc));
            return ptr;
        }

        case DiagonalizationStrategy::Subspace:
        {
            std::unique_ptr<DiagonalizationEngine> ptr(new SubspaceDiagonalizer(sc));
            return ptr;
        }

        case DiagonalizationStrategy::FirstOrder:
        {
            std::unique_ptr<DiagonalizationEngine> ptr(new FirstOrderResponseEngine(sc));
            return ptr;
        }

        default:
        {
            std::unique_ptr<DiagonalizationEngine> ptr(new FullSpaceDiagonalizer(sc));
            return ptr;
        }
    }
}

DiagonalizationStrategy strategy_from_string(const std::string& s)
{
    if (s == "full" || s == "fullspace" || s == "full_space")
        return DiagonalizationStrategy::FullSpace;
    if (s == "subspace" || s == "subspace_diag")
        return DiagonalizationStrategy::Subspace;
    if (s == "first_order" || s == "firstorder" || s == "linear_response")
        return DiagonalizationStrategy::FirstOrder;
    return DiagonalizationStrategy::FullSpace; // default fallback
}

std::string strategy_to_string(DiagonalizationStrategy s)
{
    switch (s)
    {
        case DiagonalizationStrategy::FullSpace:
            return "FullSpace";
        case DiagonalizationStrategy::Subspace:
            return "Subspace";
        case DiagonalizationStrategy::FirstOrder:
            return "FirstOrder";
        default:
            return "Unknown";
    }
}

} // namespace spinconstrain
