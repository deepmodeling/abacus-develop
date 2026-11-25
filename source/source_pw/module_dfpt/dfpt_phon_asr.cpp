// ============================================================
// Acoustic sum rule family of DFPT_Phon, the third companion
// split of dfpt_phon.cpp (ewald / elec / asr): the per-direction physical
// residual audit and the simple/crystal enforcement with the
// mass-weighted symmetric projection and the Born-charge sum
// rule (QE set_asr convention). Both act at Gamma only; the
// companion tests live in test_serial/dfpt_phon_serial_test.cpp.
// ============================================================

#include "dfpt_phon.h"
#include "dfpt_pw_data.h"

#include "source_base/timer.h"
#include "source_base/tool_title.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

namespace ModuleDFPT
{

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

} // namespace ModuleDFPT
