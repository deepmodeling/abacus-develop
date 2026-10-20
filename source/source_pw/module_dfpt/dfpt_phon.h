#ifndef DFPT_PHON_H
#define DFPT_PHON_H

#include "dfpt_pw_data.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"

#include <string>

namespace ModulePW
{
class PW_Basis;
}

namespace ModuleDFPT
{

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
class DFPT_Phon
{
  public:
    DFPT_Phon();
    ~DFPT_Phon();

    void init(UnitCell& ucell, ModulePW::PW_Basis* pw_rho, DFPT_Pert* pert);

    void assemble(int q_idx, DFPT_PW_Data& data);

    /// Fill the D[b][*] row of the electronic dynamical-matrix contribution
    /// for the converged displacement (atom_idx, dir); psi/wg are the
    /// ground-state wavefunctions and occupations. Requires a wired
    /// DFPT_Pert (init); a null pert leaves the row untouched.
    void accumulate_electron(int q_idx,
                             int atom_idx,
                             int dir,
                             const psi::Psi<std::complex<double>>& psi,
                             const ModuleBase::matrix& wg,
                             DFPT_PW_Data& data);

    void diagonalize(int q_idx, DFPT_PW_Data& data);

    /// Diagonalize the LO-TO corrected Gamma dynamical matrix (after
    /// add_loto has merged the non-analytic term into data.dynmat(0)) and
    /// store the signed frequencies (cm^-1) into data.phon_freq_loto; the
    /// uncorrected data.phon_freq(0) of the plain diagonalize stays intact.
    void diagonalize_loto(DFPT_PW_Data& data);

    /// Human-readable per-q frequency report: a header carrying the q index
    /// and the direct q coordinates plus one signed-frequency line per mode.
    /// Deterministic fixed-precision formatting; consumed by the esolver
    /// post-processing output and pinned by the format regression test.
    std::string format_q_report(int q_idx, const DFPT_PW_Data& data) const;

    /// Report of the LO-TO corrected Gamma frequencies along
    /// data.loto_dir(); returns an empty string unless the corrected
    /// frequencies have been computed (add_loto + diagonalize_loto).
    std::string format_loto_report(const DFPT_PW_Data& data) const;

    /// Non-analytic (LO-TO) term along the q->0 direction qhat (unit vector,
    /// Cartesian): D_NAC = (4 pi e^2/Omega) (qhat Z*_a)(qhat Z*_b) /
    /// (qhat eps_inf qhat) / sqrt(M_a M_b). Uses the dielectric tensor and
    /// Born charges stored in data (set by DFPT_Q0, C6).
    void add_loto(const ModuleBase::Vector3<double>& qhat, DFPT_PW_Data& data);

    /// Acoustic sum rule check at q=Gamma: max_a |sum_b D_ab| relative to
    /// the largest matrix element; returns true when below 1e-6 (or away
    /// from Gamma, where the rule does not apply).
    bool check_sum_rule(int q_idx, DFPT_PW_Data& data) const;

    /// translational sum-rule audit of the Gamma dynamical matrix
    /// (per-atom 3-vector row sums, relative to the matrix scale);
    /// diagnostics for the acoustic-branch investigation: the raw matrix
    /// and the LO-TO-corrected one are audited separately because the
    /// nonanalytic term carries its own Z* sum rule
    void audit_asr(const std::string& label, int q_idx, const DFPT_PW_Data& data) const;

    /// enforce the translational acoustic sum rule on the Gamma dynamical
    /// matrix (and, in every enforcing mode, the Z* sum rule
    /// sum_a Z*_a = 0 on the stored Born charges, like QE dynmat's
    /// set_asr). Both modes operate on the mass-reduced matrix D =
    /// Phi/sqrt(m_a m_b); the physical rule reads sum_b sqrt(m_b)
    /// D(3a+alpha, 3b+beta) = 0 for each translation direction beta
    /// independently.
    ///
    /// Both modes finish with the exact Hermitian projection D ->
    /// (1-P) D (1-P), P = sum_i w_i w_i^T, with the orthonormal physical
    /// translations w_i(3a+alpha) = sqrt(m_a) delta(alpha,i) /
    /// sqrt(sum m): the final matrix satisfies all three mass-weighted
    /// translation constraints exactly and is Hermitian.
    ///
    /// "simple": the QE legacy flavor -- the residual is first dumped
    /// onto the onsite 3x3 blocks, D(a,a)[alpha][beta] =
    /// -sum_{b!=a} sqrt(m_b/m_a) D(a,alpha; b,beta), then the projection
    /// absorbs whatever onsite-only correction could not (a purely
    /// onsite Hermitian fix does not exist when the residual has a
    /// cross-Cartesian antisymmetric part).
    ///
    /// "crystal": the optimized minimal-norm fix -- the projection alone,
    /// which also keeps symmetry-equivalent atoms equivalent.
    ///
    /// Away from Gamma the call is a no-op.
    void apply_asr(const std::string& asr_type, int q_idx, DFPT_PW_Data& data);

    /// Ewald ion-ion force constants C^ewald_ab(q) (G-space + real-space +
    /// Gaussian self term), mass-reduced by 1/sqrt(M_a M_b). A stateless
    /// building block exposed publicly so the serial analytic tests can
    /// validate it directly (no internal state is touched).
    void ion_ion(const ModuleBase::Vector3<double>& q_frac, ModuleBase::ComplexMatrix& dyn);

    /// read-only view of the accumulated electronic dynamical-matrix rows
    /// (see dynmat_accum_); consumed by the serial analytic tests.
    const ModuleBase::ComplexMatrix& dynmat_accum() const
    {
        return dynmat_accum_;
    }

  private:
    UnitCell* ucell_ = nullptr;
    ModulePW::PW_Basis* pw_rho_ = nullptr;
    DFPT_Pert* pert_ = nullptr;

    double ewald_alpha_ = 0.0;
    double ewald_rcut_ = 0.0;

    /// same-atom anharmonic term <psi | d2_ab V_ext | psi> accumulated into
    /// the (rowb, cola) entries of dynmat_accum_ (upper triangle only; the
    /// Hermitian partner is added here from conj(d2sum))
    void accum_d2_same_atom(int atom_idx,
                            int dir,
                            int iat,
                            int idir,
                            int rowb,
                            const psi::Psi<std::complex<double>>& psi,
                            const ModuleBase::matrix& wg);

    /// DFT+U contribution to the dynamical matrix (U0 reservation).
    void dftu_onsite(int q_idx, DFPT_PW_Data& data);

    /// accumulated electronic rows of the current q (merged and cleared by
    /// assemble)
    ModuleBase::ComplexMatrix dynmat_accum_;
    int accum_q_ = -1;
};

} // namespace ModuleDFPT

#endif // DFPT_PHON_H
