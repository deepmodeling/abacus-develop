// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_PHON_H
#define DFPT_PHON_H

#include "dfpt_pw_data.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"

namespace ModulePW {
class PW_Basis;
}

namespace ModuleDFPT {

class DFPT_Pert;

/**
 * @brief Dynamical matrix of DFPT (C5).
 *
 * assemble() merges the Ewald ion-ion force constants with the electronic
 * contribution accumulated by accumulate_electron(); diagonalize() solves
 * the complex Hermitian eigenproblem through LapackConnector::zheev and
 * stores frequencies as signed values omega_i = sgn(e_i) sqrt(|e_i|) in
 * cm^-1 (negative = imaginary frequency).
 *
 * Electronic contribution (2n+1 theorem, insulating case):
 *   D_ab = D^Ewald_ab
 *        + 2 sum_kn wg Re <dpsi^b_nk | dV^a_ext | psi_nk>
 *        + sum_kn wg <psi_nk | d2_ab V_ext | psi_nk>   (same atom a,b only)
 * with dpsi^b the converged Sternheimer solution for displacement b and
 * dV^a_ext the BARE first-order external potential of displacement a; the
 * row D[b][*] is filled right after displacement b converges, so the dpsi
 * storage never needs a direction dimension (data-layer refactor reserved
 * for phase B).
 */
class DFPT_Phon {
public:
    DFPT_Phon();
    ~DFPT_Phon();
    
    void init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, DFPT_Pert* pert);
    
    void assemble(int q_idx, DFPT_PW_Data& data);

    /// Fill the D[b][*] row of the electronic dynamical-matrix contribution
    /// for the converged displacement (atom_idx, dir); psi/wg are the
    /// ground-state wavefunctions and occupations. Requires a wired
    /// DFPT_Pert (init); a null pert leaves the row untouched.
    void accumulate_electron(int q_idx, int atom_idx, int dir,
                             const psi::Psi<std::complex<double>>& psi,
                             const ModuleBase::matrix& wg, DFPT_PW_Data& data);
    
    void diagonalize(int q_idx, DFPT_PW_Data& data);
    
    /// Non-analytic (LO-TO) term along the q->0 direction qhat (unit vector,
    /// Cartesian): D_NAC = (4 pi e^2/Omega) (qhat Z*_a)(qhat Z*_b) /
    /// (qhat eps_inf qhat) / sqrt(M_a M_b). Uses the dielectric tensor and
    /// Born charges stored in data (set by DFPT_Q0, C6).
    void add_loto(const ModuleBase::Vector3<double>& qhat, DFPT_PW_Data& data);
    
    /// Acoustic sum rule check at q=Gamma: max_a |sum_b D_ab| relative to
    /// the largest matrix element; returns true when below 1e-6 (or away
    /// from Gamma, where the rule does not apply).
    bool check_sum_rule(int q_idx, DFPT_PW_Data& data) const;

private:
    UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    DFPT_Pert* pert_ = nullptr;
    
    double ewald_alpha_ = 0.0;
    double ewald_rcut_ = 0.0;
    
    /// Ewald ion-ion force constants C^ewald_ab(q) (G-space + real-space +
    /// Gaussian self term), mass-reduced by 1/sqrt(M_a M_b).
    void ion_ion(const ModuleBase::Vector3<double>& q_frac, ModuleBase::ComplexMatrix& dyn);
    
    /// DFT+U contribution to the dynamical matrix (U0 reservation).
    void dftu_onsite(int q_idx, DFPT_PW_Data& data);

    /// accumulated electronic rows of the current q (merged and cleared by
    /// assemble)
    ModuleBase::ComplexMatrix dynmat_accum_;
    int accum_q_ = -1;
};

} // namespace ModuleDFPT

#endif // DFPT_PHON_H
