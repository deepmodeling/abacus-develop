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
