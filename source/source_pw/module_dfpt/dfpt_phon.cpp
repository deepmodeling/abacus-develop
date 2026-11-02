// The Ewald ion-ion part of ion_ion lives in dfpt_phon_ewald.cpp and
// the electronic (2n+1) part of accumulate_electron in
// dfpt_phon_elec.cpp.

#include "dfpt_phon.h"

#include "source_base/constants.h"
#include "source_base/global_function.h"
#include "source_base/module_external/lapack_connector.h"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "source_base/timer.h"
#include "source_base/tool_title.h"

namespace ModuleDFPT
{

DFPT_Phon::DFPT_Phon()
{
}

DFPT_Phon::~DFPT_Phon()
{
}

namespace
{

// signed frequencies: omega = sgn(e) sqrt(|e|), converted to cm^-1
// sqrt(Ry/(bohr^2 amu)) in cm^-1 = sqrt(RYDBERG_SI/amu_kg)/(bohr*2pi*c)
std::vector<double> signed_freqs_cm1(const std::vector<double>& eigs)
{
    ModuleBase::TITLE("DFPT_Phon", "signed_freqs_cm1");
    ModuleBase::timer::start("DFPT_Phon", "signed_freqs_cm1");
    const double amu_kg = 1.0e-3 / ModuleBase::NA;
    const double light_speed_cgs = 2.99792458e10; // cm/s, exact SI value
    const double ry_bohr2_amu_to_cm1 = std::sqrt(ModuleBase::RYDBERG_SI / amu_kg)
                                       / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::TWO_PI * light_speed_cgs);
    std::vector<double> freq(eigs.size(), 0.0);
    for (size_t i = 0; i < eigs.size(); ++i)
    {
        freq[i] = ((eigs[i] >= 0.0) ? 1.0 : -1.0) * std::sqrt(std::abs(eigs[i])) * ry_bohr2_amu_to_cm1;
    }
    ModuleBase::timer::end("DFPT_Phon", "signed_freqs_cm1");
    return freq;
}

} // namespace

void DFPT_Phon::init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, DFPT_Pert* pert)
{
    ModuleBase::TITLE("DFPT_Phon", "init");
    ModuleBase::timer::start("DFPT_Phon", "init");
    ucell_ = &ucell;
    pw_rho_ = pw_rho;
    pert_ = pert;
    ModuleBase::timer::end("DFPT_Phon", "init");
}

void DFPT_Phon::set_rho(const DFPT_Rho* rho)
{
    rho_ = rho;
}

void DFPT_Phon::set_vxc_r(const std::vector<double>* vxc_r)
{
    vxc_r_ = vxc_r;
}

void DFPT_Phon::set_dvxc_core_table(const std::vector<std::vector<std::complex<double>>>* table)
{
    dvxc_core_ = table;
}

// ---------------------------------------------------------------------------
// assemble / diagonalize / LO-TO / sum rule
// ---------------------------------------------------------------------------

void DFPT_Phon::assemble(int q_idx, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Phon", "assemble");
    ModuleBase::timer::start("DFPT_Phon", "assemble");
    if (ucell_ == nullptr)
    {
        ModuleBase::timer::end("DFPT_Phon", "assemble");
        return;
    }
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    ModuleBase::ComplexMatrix dyn(nat3, nat3, true);
    if (pw_rho_ != nullptr)
    {
        ion_ion(data.get_qvec(q_idx), dyn);
    }
    if (accum_q_ == q_idx && dynmat_accum_.nr == nat3)
    {
        for (int i = 0; i < nat3; ++i)
        {
            for (int j = 0; j < nat3; ++j)
            {
                dyn(i, j) += dynmat_accum_(i, j);
            }
        }
    }
    // DFT+U dynamical-matrix term (U0 reservation, implemented with the C7/U1
    // Plus_U wiring): sum_nk w_nk [<psi|dV_U|dpsi> + frozen second-order term].
    if (data.with_u())
    {
        dftu_onsite(q_idx, data);
    }
    // Hermitian symmetrization (rows filled by independent solves)
    for (int i = 0; i < nat3; ++i)
    {
        for (int j = i + 1; j < nat3; ++j)
        {
            const std::complex<double> avg = 0.5 * (dyn(i, j) + std::conj(dyn(j, i)));
            dyn(i, j) = avg;
            dyn(j, i) = std::conj(avg);
        }
    }
    data.set_dynmat(q_idx, dyn);
    // stash the raw pre-ASR pre-LO-TO copy for the dyn file writer;
    // later set_dynmat calls (apply_asr, add_loto) update the working
    // matrix but must not overwrite this snapshot
    data.set_dynmat_raw(q_idx, dyn);
    dynmat_accum_ = ModuleBase::ComplexMatrix();
    accum_q_ = -1;
    ModuleBase::timer::end("DFPT_Phon", "assemble");
}

void DFPT_Phon::diagonalize(int q_idx, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Phon", "diagonalize");
    ModuleBase::timer::start("DFPT_Phon", "diagonalize");
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    ModuleBase::ComplexMatrix dyn = data.get_dynmat(q_idx);
    if (dyn.nr != nat3)
    {
        ModuleBase::timer::end("DFPT_Phon", "diagonalize");
        return;
    }

    // eigenvalues of the complex Hermitian dynamical matrix (Ry/bohr^2/amu)
    std::vector<double> w(nat3, 0.0);
    std::vector<double> rwork(std::max(1, 3 * nat3 - 2), 0.0);
    std::vector<std::complex<double>> work(1);
    int info = 0;
    LapackConnector::zheev('N', 'U', nat3, dyn, nat3, w.data(), work.data(), -1, rwork.data(), &info);
    work.resize(std::max(1, static_cast<int>(work[0].real())));
    LapackConnector::zheev('N',
                           'U',
                           nat3,
                           dyn,
                           nat3,
                           w.data(),
                           work.data(),
                           static_cast<int>(work.size()),
                           rwork.data(),
                           &info);

    // signed frequencies: omega = sgn(e) sqrt(|e|), converted to cm^-1
    // sqrt(Ry/(bohr^2 amu)) in cm^-1 = sqrt(RYDBERG_SI/amu_kg)/(bohr*2pi*c)
    const double amu_kg = 1.0e-3 / ModuleBase::NA;
    const double light_speed_cgs = 2.99792458e10; // cm/s, exact SI value
    const double ry_bohr2_amu_to_cm1 = std::sqrt(ModuleBase::RYDBERG_SI / amu_kg)
                                       / (ModuleBase::BOHR_RADIUS_SI * ModuleBase::TWO_PI * light_speed_cgs);
    std::vector<double> freq(nat3, 0.0);
    for (int i = 0; i < nat3; ++i)
    {
        const double e = w[i];
        freq[i] = ((e >= 0.0) ? 1.0 : -1.0) * std::sqrt(std::abs(e)) * ry_bohr2_amu_to_cm1;
    }
    data.set_phon_freq(q_idx, freq);
    ModuleBase::timer::end("DFPT_Phon", "diagonalize");
}

void DFPT_Phon::audit_asr(const std::string& label, int q_idx, const DFPT_PW_Data& data) const
{
    ModuleBase::TITLE("DFPT_Phon", "audit_asr");
    ModuleBase::timer::start("DFPT_Phon", "audit_asr");
    const double gamma_tol = 1.0e-8; ///< empirical parameter: fractional-q Gamma tolerance
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    const bool at_gamma = std::abs(q_frac.x) <= gamma_tol && std::abs(q_frac.y) <= gamma_tol
                          && std::abs(q_frac.z) <= gamma_tol;
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    const ModuleBase::ComplexMatrix dyn = data.get_dynmat(q_idx);
    if (!at_gamma || dyn.nr != nat3)
    {
        ModuleBase::timer::end("DFPT_Phon", "audit_asr");
        return;
    }
    double max_elem = 0.0;
    for (int i = 0; i < nat3; ++i)
    {
        for (int j = 0; j < nat3; ++j)
        {
            max_elem = std::max(max_elem, std::abs(dyn(i, j)));
        }
    }
    // per-atom PHYSICAL translational row sums: the stored matrix is
    // mass-reduced (D = Phi/sqrt(m_a m_b)), so the Phi row sums read
    // s(ia)[da] = sum_jb sqrt(m_jb) D(3ia+da, jb) -- the same quantity
    // QE's unreduced dyn rows sum to
    std::vector<double> sqm(nat, 0.0);
    for (int ia = 0; ia < nat; ++ia)
    {
        sqm[ia] = std::sqrt(ucell_->atoms[ucell_->iat2it[ia]].mass);
    }
    std::cout << "DFPT ASR audit (" << label << "):";
    double worst = 0.0;
    int worst_atom = -1;
    for (int ia = 0; ia < nat; ++ia)
    {
        // per-DIRECTION residuals R(a,alpha,i) = sum_b sqrt(m_b) D(3a+alpha,
        // 3b+i): summing over atoms at a FIXED translation direction i
        // before combining norms, so opposite x/y residuals cannot cancel
        double s2 = 0.0;
        for (int da = 0; da < 3; ++da)
        {
            for (int db = 0; db < 3; ++db)
            {
                std::complex<double> s(0.0, 0.0);
                for (int ib = 0; ib < nat; ++ib)
                {
                    s += sqm[ib] * dyn(3 * ia + da, 3 * ib + db);
                }
                s2 += std::norm(s);
            }
        }
        const double snorm = std::sqrt(s2);
        std::cout << " atom" << ia + 1 << "=" << std::setprecision(3) << std::scientific << snorm;
        if (snorm > worst)
        {
            worst = snorm;
            worst_atom = ia;
        }
    }
    std::cout << std::defaultfloat << std::setprecision(6);
    if (max_elem > 0.0 && worst_atom >= 0)
    {
        std::cout << " | worst rel = " << std::setprecision(3) << std::scientific << worst / max_elem
                  << std::defaultfloat << std::setprecision(6);
    }
    std::cout << std::endl;
    ModuleBase::timer::end("DFPT_Phon", "audit_asr");
}

void DFPT_Phon::apply_asr(const std::string& asr_type, int q_idx, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Phon", "apply_asr");
    ModuleBase::timer::start("DFPT_Phon", "apply_asr");
    const double gamma_tol = 1.0e-8; ///< empirical parameter: fractional-q Gamma tolerance
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    const bool at_gamma = std::abs(q_frac.x) <= gamma_tol && std::abs(q_frac.y) <= gamma_tol
                          && std::abs(q_frac.z) <= gamma_tol;
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    ModuleBase::ComplexMatrix dyn = data.get_dynmat(q_idx);
    if (!at_gamma || dyn.nr != nat3)
    {
        std::cout << "DFPT ASR (" << asr_type << ") skipped: not a Gamma matrix" << std::endl;
        ModuleBase::timer::end("DFPT_Phon", "apply_asr");
        return;
    }

    // The stored matrix is mass-reduced, D = Phi / sqrt(m_a m_b) (amu masses
    // from the cell); the PHYSICAL sum rule acts on Phi rows,
    // sum_b sqrt(m_b) D(3a+da, 3b+db) = 0, and the acoustic eigenvectors of
    // D are v(3a+da) = sqrt(m_a) e_da (uniform physical displacement). QE
    // dynmat works on the unreduced Phi, so its unweighted rows are this
    // same rule; here the constraint vectors carry the sqrt(m) weights.
    std::vector<double> sqm(nat, 0.0);
    double mtot = 0.0;
    for (int ia = 0; ia < nat; ++ia)
    {
        sqm[ia] = std::sqrt(ucell_->atoms[ucell_->iat2it[ia]].mass);
        mtot += ucell_->atoms[ucell_->iat2it[ia]].mass;
    }
    // per-DIRECTION physical residual R(a,alpha,i) = sum_b sqrt(m_b)
    // D(3a+alpha, 3b+i): summing over atoms at a fixed translation
    // direction i so opposite x/y residuals cannot cancel
    const auto physical_residual_norm = [&](int ia) {
        double s2 = 0.0;
        for (int da = 0; da < 3; ++da)
        {
            for (int db = 0; db < 3; ++db)
            {
                std::complex<double> s(0.0, 0.0);
                for (int ib = 0; ib < nat; ++ib)
                {
                    s += sqm[ib] * dyn(3 * ia + da, 3 * ib + db);
                }
                s2 += std::norm(s);
            }
        }
        return std::sqrt(s2);
    };

    // violation BEFORE the correction (physical residuals)
    double worst_before = 0.0;
    for (int ia = 0; ia < nat; ++ia)
    {
        worst_before = std::max(worst_before, physical_residual_norm(ia));
    }

    // the crystal projection D -> (1-P) D (1-P) with the orthonormal
    // physical translations w_i(3a+da) = sqrt(m_a) delta(da,i)/sqrt(sum m):
    // zeroes every per-direction physical residual exactly, keeps
    // Hermiticity and symmetry-equivalent atoms equivalent, and leaves the
    // acoustic eigenvectors in the null space so the LO-TO term (which
    // decouples from them through sum_a Z*_a = 0) cannot lift the acoustic
    // frequencies. Both 'simple' and 'crystal' finish with this projection:
    // 'simple' first dumps the residual onto the onsite blocks (the QE
    // legacy flavor), then the projection restores Hermiticity and removes
    // whatever onsite-only correction could not absorb (cross-Cartesian
    // antisymmetric residuals have no exact Hermitian onsite-only fix)
    const auto apply_crystal_projection = [&]() {
        std::vector<std::complex<double>> s(3 * nat3, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> t(3 * nat3, std::complex<double>(0.0, 0.0));
        ModuleBase::ComplexMatrix u(3, 3, true);
        for (int na = 0; na < nat; ++na)
        {
            for (int da = 0; da < 3; ++da)
            {
                for (int nb = 0; nb < nat; ++nb)
                {
                    for (int db = 0; db < 3; ++db)
                    {
                        const std::complex<double> d = dyn(3 * na + da, 3 * nb + db);
                        s[da * nat3 + 3 * nb + db] += sqm[na] * d;
                        t[(3 * na + da) * 3 + db] += sqm[nb] * d;
                        u(da, db) += sqm[na] * sqm[nb] * d;
                    }
                }
            }
        }
        for (int na = 0; na < nat; ++na)
        {
            for (int da = 0; da < 3; ++da)
            {
                for (int nb = 0; nb < nat; ++nb)
                {
                    for (int db = 0; db < 3; ++db)
                    {
                        std::complex<double> d = dyn(3 * na + da, 3 * nb + db);
                        d -= sqm[na] / mtot * s[da * nat3 + 3 * nb + db];
                        d -= sqm[nb] / mtot * t[(3 * na + da) * 3 + db];
                        d += sqm[na] * sqm[nb] / (mtot * mtot) * u(da, db);
                        dyn(3 * na + da, 3 * nb + db) = d;
                    }
                }
            }
        }
    };

    if (asr_type == "simple")
    {
        // legacy QE 'simple' translated to the reduced matrix: dump the
        // per-direction residual onto the onsite 3x3 block,
        // D(a,a)[da][db] = -sum_{b!=a} sqrt(m_b/m_a) D(a,da; b,db);
        // the onsite dump may break Hermiticity (cross-Cartesian parts),
        // so symmetrize, then the crystal projection below removes the
        // residual the symmetrization reintroduced (at most half the
        // anti-Hermitian part) and guarantees the final matrix is
        // Hermitian and exactly satisfies all three translation constraints
        for (int na = 0; na < nat; ++na)
        {
            for (int da = 0; da < 3; ++da)
            {
                for (int db = 0; db < 3; ++db)
                {
                    std::complex<double> s(0.0, 0.0);
                    for (int nb = 0; nb < nat; ++nb)
                    {
                        if (nb != na)
                        {
                            s += sqm[nb] * dyn(3 * na + da, 3 * nb + db);
                        }
                    }
                    dyn(3 * na + da, 3 * na + db) = -s / sqm[na];
                }
            }
        }
        // Hermitian re-symmetrization of the onsite dump (the projection
        // is Hermitian-preserving only on Hermitian input); the j > i loop
        // alone leaves the diagonal untouched, and the onsite dump may
        // have written complex values there -- take the real part so the
        // matrix is genuinely Hermitian before the projection
        for (int i = 0; i < nat3; ++i)
        {
            for (int j = i + 1; j < nat3; ++j)
            {
                const std::complex<double> avg = 0.5 * (dyn(i, j) + std::conj(dyn(j, i)));
                dyn(i, j) = avg;
                dyn(j, i) = std::conj(avg);
            }
            dyn(i, i) = std::complex<double>(dyn(i, i).real(), 0.0);
        }
        apply_crystal_projection();
    }
    else if (asr_type == "crystal")
    {
        apply_crystal_projection();
    }
    else
    {
        std::cout << "DFPT ASR: unknown type '" << asr_type << "', nothing applied" << std::endl;
        ModuleBase::timer::end("DFPT_Phon", "apply_asr");
        return;
    }

    // Hermitian re-symmetrization (rows are corrected independently)
    for (int i = 0; i < nat3; ++i)
    {
        for (int j = i + 1; j < nat3; ++j)
        {
            const std::complex<double> avg = 0.5 * (dyn(i, j) + std::conj(dyn(j, i)));
            dyn(i, j) = avg;
            dyn(j, i) = std::conj(avg);
        }
    }

    // Z* sum rule (both enforcing modes, QE set_asr convention): shift by
    // the per-component mean over atoms, sum_a Z*_a -> 0
    double zworst_before = 0.0;
    ModuleBase::ComplexMatrix zsum(3, 3, true);
    int n_born = 0;
    for (int ia = 0; ia < nat; ++ia)
    {
        const ModuleBase::matrix za = data.get_born(ia);
        if (za.nr != 3 || za.nc != 3)
        {
            continue;
        }
        ++n_born;
        for (int a = 0; a < 3; ++a)
        {
            for (int b = 0; b < 3; ++b)
            {
                zsum(a, b) += std::complex<double>(za(a, b), 0.0);
            }
        }
    }
    if (n_born == nat)
    {
        for (int a = 0; a < 3; ++a)
        {
            for (int b = 0; b < 3; ++b)
            {
                zworst_before = std::max(zworst_before, std::abs(zsum(a, b)));
            }
        }
        for (int ia = 0; ia < nat; ++ia)
        {
            ModuleBase::matrix za = data.get_born(ia);
            for (int a = 0; a < 3; ++a)
            {
                for (int b = 0; b < 3; ++b)
                {
                    za(a, b) -= zsum(a, b).real() / nat;
                }
            }
            data.set_born(ia, za);
        }
    }

    data.set_dynmat(q_idx, dyn);

    // violation AFTER the correction (physical residuals)
    double worst_after = 0.0;
    for (int ia = 0; ia < nat; ++ia)
    {
        worst_after = std::max(worst_after, physical_residual_norm(ia));
    }
    std::cout << "DFPT ASR (" << asr_type << ") applied: worst physical residual " << std::setprecision(3)
              << std::scientific << worst_before << " -> " << worst_after;
    if (n_born == nat)
    {
        std::cout << " | max |sum_a Z*_a| " << zworst_before << " -> 0";
    }
    std::cout << std::defaultfloat << std::setprecision(6) << std::endl;
    ModuleBase::timer::end("DFPT_Phon", "apply_asr");
}

void DFPT_Phon::add_loto(const ModuleBase::Vector3<double>& qhat, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Phon", "add_loto");
    ModuleBase::timer::start("DFPT_Phon", "add_loto");
    const int nat = ucell_->nat;
    const int nat3 = 3 * nat;
    ModuleBase::ComplexMatrix dyn = data.get_dynmat(0);
    const double qeq_tol = 1.0e-10; ///< empirical parameter: |q eps q| floor for the non-polar-direction skip
    if (dyn.nr != nat3)
    {
        ModuleBase::timer::end("DFPT_Phon", "add_loto");
        return;
    }
    const ModuleBase::matrix eps = data.get_dielectric();
    if (eps.nr != 3 || eps.nc != 3)
    {
        ModuleBase::timer::end("DFPT_Phon", "add_loto");
        return; // no dielectric tensor stored yet (C6 not run)
    }
    const double qeq = qhat.x * (qhat.x * eps(0, 0) + qhat.y * eps(1, 0) + qhat.z * eps(2, 0))
                       + qhat.y * (qhat.x * eps(0, 1) + qhat.y * eps(1, 1) + qhat.z * eps(2, 1))
                       + qhat.z * (qhat.x * eps(0, 2) + qhat.y * eps(1, 2) + qhat.z * eps(2, 2));
    if (std::abs(qeq) < qeq_tol)
    {
        ModuleBase::timer::end("DFPT_Phon", "add_loto");
        return;
    }
    const double pref = ModuleBase::FOUR_PI * ModuleBase::e2 / ucell_->omega / qeq;
    for (int ia = 0; ia < nat; ++ia)
    {
        const double ma = ucell_->atoms[ucell_->iat2it[ia]].mass;
        const ModuleBase::matrix za = data.get_born(ia);
        if (za.nr != 3 || za.nc != 3)
        {
            continue;
        }
        for (int ib = 0; ib < nat; ++ib)
        {
            const double mb = ucell_->atoms[ucell_->iat2it[ib]].mass;
            const ModuleBase::matrix zb = data.get_born(ib);
            for (int da = 0; da < 3; ++da)
            {
                // (qhat Z*_a)_da = sum_gamma qhat_gamma Z_a(da,gamma)
                const double qza = qhat.x * za(da, 0) + qhat.y * za(da, 1) + qhat.z * za(da, 2);
                for (int db = 0; db < 3; ++db)
                {
                    const double qzb = qhat.x * zb(db, 0) + qhat.y * zb(db, 1) + qhat.z * zb(db, 2);
                    dyn(3 * ia + da, 3 * ib + db) += pref * qza * qzb / std::sqrt(ma * mb);
                }
            }
        }
    }
    data.set_dynmat(0, dyn);
    ModuleBase::timer::end("DFPT_Phon", "add_loto");
}

void DFPT_Phon::diagonalize_loto(DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Phon", "diagonalize_loto");
    ModuleBase::timer::start("DFPT_Phon", "diagonalize_loto");
    const int nat3 = 3 * ucell_->nat;
    // the stored Gamma matrix already carries the non-analytic term added
    // by add_loto; the copy below is destroyed by the solver, the stored
    // one stays available for the plain report
    ModuleBase::ComplexMatrix dyn = data.get_dynmat(0);
    if (dyn.nr != nat3)
    {
        ModuleBase::timer::end("DFPT_Phon", "diagonalize_loto");
        return;
    }
    std::vector<double> w(nat3, 0.0);
    std::vector<double> rwork(std::max(1, 3 * nat3 - 2), 0.0);
    std::vector<std::complex<double>> work(1);
    int info = 0;
    LapackConnector::zheev('N', 'U', nat3, dyn, nat3, w.data(), work.data(), -1, rwork.data(), &info);
    work.resize(std::max(1, static_cast<int>(work[0].real())));
    LapackConnector::zheev('N',
                           'U',
                           nat3,
                           dyn,
                           nat3,
                           w.data(),
                           work.data(),
                           static_cast<int>(work.size()),
                           rwork.data(),
                           &info);
    data.set_phon_freq_loto(signed_freqs_cm1(w));
    ModuleBase::timer::end("DFPT_Phon", "diagonalize_loto");
}

std::string DFPT_Phon::format_q_report(int q_idx, const DFPT_PW_Data& data) const
{
    ModuleBase::TITLE("DFPT_Phon", "format_q_report");
    ModuleBase::timer::start("DFPT_Phon", "format_q_report");
    const ModuleBase::Vector3<double> qd = data.get_qvec(q_idx);
    const std::vector<double> freq = data.get_phon_freq(q_idx);
    std::ostringstream os;
    os << " DFPT phonon frequencies at q #" << q_idx << " = (" << std::fixed << std::setprecision(6) << qd.x << " "
       << qd.y << " " << qd.z << ") (direct) in cm^-1:" << "\n";
    for (size_t im = 0; im < freq.size(); ++im)
    {
        os << "   mode " << std::setw(3) << im << " : " << std::fixed << std::setprecision(6) << freq[im] << " cm^-1"
           << "\n";
    }
    ModuleBase::timer::end("DFPT_Phon", "format_q_report");
    return os.str();
}

std::string DFPT_Phon::format_loto_report(const DFPT_PW_Data& data) const
{
    ModuleBase::TITLE("DFPT_Phon", "format_loto_report");
    ModuleBase::timer::start("DFPT_Phon", "format_loto_report");
    const std::vector<double> freq = data.get_phon_freq_loto();
    if (freq.empty())
    {
        ModuleBase::timer::end("DFPT_Phon", "format_loto_report");
        return std::string();
    }
    const ModuleBase::Vector3<double> dir = data.get_loto_dir();
    std::ostringstream os;
    os << " DFPT LO-TO corrected frequencies at q #0 along q->0 direction (" << std::fixed << std::setprecision(6)
       << dir.x << " " << dir.y << " " << dir.z << ") in cm^-1:" << "\n";
    for (size_t im = 0; im < freq.size(); ++im)
    {
        os << "   mode " << std::setw(3) << im << " : " << std::fixed << std::setprecision(6) << freq[im] << " cm^-1"
           << "\n";
    }
    ModuleBase::timer::end("DFPT_Phon", "format_loto_report");
    return os.str();
}

bool DFPT_Phon::check_sum_rule(int q_idx, DFPT_PW_Data& data) const
{
    ModuleBase::TITLE("DFPT_Phon", "check_sum_rule");
    ModuleBase::timer::start("DFPT_Phon", "check_sum_rule");
    const double gamma_tol = 1.0e-8;       ///< empirical parameter: fractional-q Gamma tolerance
    const double dyn_zero_floor = 1.0e-12; ///< empirical parameter: dynamical matrix treated as zero
    const double asr_rel_tol = 1.0e-6;     ///< empirical parameter: tolerated relative column-sum error
    const ModuleBase::Vector3<double> q_frac = data.get_qvec(q_idx);
    if (std::abs(q_frac.x) > gamma_tol || std::abs(q_frac.y) > gamma_tol || std::abs(q_frac.z) > gamma_tol)
    {
        ModuleBase::timer::end("DFPT_Phon", "check_sum_rule");
        return true; // only applies at Gamma
    }
    const int nat3 = 3 * ucell_->nat;
    ModuleBase::ComplexMatrix dyn = data.get_dynmat(q_idx);
    if (dyn.nr != nat3)
    {
        ModuleBase::timer::end("DFPT_Phon", "check_sum_rule");
        return false;
    }
    double max_elem = 0.0;
    for (int i = 0; i < nat3; ++i)
    {
        for (int j = 0; j < nat3; ++j)
        {
            max_elem = std::max(max_elem, std::abs(dyn(i, j)));
        }
    }
    if (max_elem < dyn_zero_floor)
    {
        ModuleBase::timer::end("DFPT_Phon", "check_sum_rule");
        return true;
    }
    for (int i = 0; i < nat3; ++i)
    {
        std::complex<double> colsum(0.0, 0.0);
        for (int j = 0; j < nat3; ++j)
        {
            colsum += dyn(i, j);
        }
        if (std::abs(colsum) > asr_rel_tol * max_elem)
        {
            ModuleBase::timer::end("DFPT_Phon", "check_sum_rule");
            return false;
        }
    }
    ModuleBase::timer::end("DFPT_Phon", "check_sum_rule");
    return true;
}

void DFPT_Phon::dftu_onsite(int q_idx, DFPT_PW_Data& data)
{
    ModuleBase::TITLE("DFPT_Phon", "dftu_onsite");
    ModuleBase::timer::start("DFPT_Phon", "dftu_onsite");
    // Reserved DFT+U contribution to the dynamical matrix (U0).
    // The physical implementation lands with the Plus_U production wiring:
    //   sum_nk w_nk [ <psi|dV_U|dpsi> + frozen second-order U term
    //   (~ becp * V_eff * dbecp_f contractions) ], accumulated into the
    //   dynamical matrix. dV_U itself is assembled by DFPT_Pert::build_dv_u.
    (void)q_idx;
    (void)data;
}

} // namespace ModuleDFPT
