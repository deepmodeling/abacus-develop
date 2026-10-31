#ifndef DFPT_PW_H
#define DFPT_PW_H

#include "source_base/matrix.h"
#include "source_base/vector3.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"

#include <memory>
#include <string>
#include <vector>

class Plus_U_Base;
class Structure_Factor;

namespace ModulePW
{
class PW_Basis;
class PW_Basis_K;
} // namespace ModulePW

namespace ModuleDFPT
{
class DFPT_PW_Data;
}

namespace ModuleDFPT
{

class XC_First_Order;

/// Bundled initialisation context for DFPT_PW (see dfpt_pw_impl.h for the
/// full field-level doxygen). Forward-declared here so callers can build a
/// struct aggregate without pulling the heavy impl header; the detailed
/// field docs live alongside the private impl header that actually uses
/// each field.
struct DFPT_PW_InitContext
{
    UnitCell* ucell;
    const psi::Psi<std::complex<double>>* psi;
    ModulePW::PW_Basis* pw_rho;
    ModulePW::PW_Basis_K* pw_wfc;
    Structure_Factor* sf;
    const std::vector<double>* veff_r;
    const ModuleBase::matrix* wg;
    const ModuleBase::matrix* eig;
    const XC_First_Order* xc;
    double nelec;
    double ecutwfc;
    const Plus_U_Base* dftu;
    /// MPI rank of the running process (from the runtime context, not
    /// read from any global inside module_dfpt); consumed by the dyn
    /// file writer for the rank-0 output gate. No default member
    /// initializer: C++11 aggregate rules require it.
    int myrank;
    /// ground-state XC potential on the rho grid (valence + core charge
    /// input, the same density the kernel linearizes around); consumed by
    /// the nlcc diagonal 2n+1 term (QE dynmatcc). Null disables that term.
    /// No default member initializer: C++11 aggregate rules require it.
    const std::vector<double>* vxc_r;
};

/**
 * @brief Density-functional perturbation theory driver (plane waves).
 *
 * C7 wiring: init receives the converged ground state (psi, wg, eig), the
 * shared-grid plane-wave bases, the real-space effective potential and the
 * first-order XC kernel contract; run() then drives, per irreducible q, the
 * per-displacement self-consistent Sternheimer cycle
 *   build_dv (bare external)
 *   -> [ dv_sc = v_hartree_q(drho_in) + xc_->apply(drho_in)
 *        -> Sternheimer solve of (H(k+q) - eps_n) P_c dpsi = -P_c (dV_ext
 *           + dV_sc)|psi_n> with the k+q occupied states as the projector
 *        -> compute_drho -> mix_drho ]*
 *   -> accumulate_electron -> assemble -> diagonalize (+ LO-TO at q = 0).
 * With null bases (design-phase skeleton) run() keeps the documented
 * first-iteration-converged fallback of the irrep bookkeeping loop.
 */
class DFPT_PW
{
  public:
    class Impl; // pimpl forward declaration (kept in public section so the
                // private nested class can be named as DFPT_PW::Impl from
                // outside translation units that include the impl header).

    DFPT_PW();
    ~DFPT_PW();

    /// Package-and-forward convenience wrapper: the former 12-argument
    /// signature is retained for backward compatibility with the small
    /// number of call sites (esolver_dfpt_pw.cpp + three test fixtures),
    /// and the actual validation/submodule wiring happens in the
    /// InitContext overload below. Keeping the thin wrapper inline avoids
    /// a separate TU and gives the call-site aggregate initialization the
    /// same performance as a direct call.
    void init(UnitCell& ucell,
              const psi::Psi<std::complex<double>>& psi,
              ModulePW::PW_Basis* pw_rho,
              ModulePW::PW_Basis_K* pw_wfc,
              Structure_Factor* sf,
              const std::vector<double>& veff_r,
              const ModuleBase::matrix& wg,
              const ModuleBase::matrix& eig,
              const XC_First_Order* xc,
              double nelec,
              double ecutwfc,
              const Plus_U_Base* dftu,
              int myrank)
    {
        const DFPT_PW_InitContext ctx{&ucell, &psi, pw_rho, pw_wfc, sf, &veff_r, &wg, &eig, xc, nelec, ecutwfc, dftu, myrank, nullptr};
        init(ctx);
    }

    /// Single-context init carrying the twelve parameters as named fields
    /// so the function body stays under the coding-rule parameter-count
    /// budget. Semantics are identical to the overload above.
    void init(const DFPT_PW_InitContext& ctx);

    void run();

    /// DFT+U reservation accessors (U0): with_u() reports whether a DFT+U
    /// provider is wired (dft_plus_u enabled upstream); u_active() further
    /// requires the provider to be usable (occupation matrices initialized).
    /// init() rejects a wired provider explicitly: every DFPT U hook is a
    /// no-op reservation, so a DFT+U ground state must not run DFPT until
    /// the first-order U response (U1) is implemented.
    bool get_with_u() const;
    bool get_u_active() const;

    std::vector<double> get_phonon_freq(int q_idx) const;

    ModuleBase::matrix get_dielectric_tensor() const;

    ModuleBase::matrix get_born_charges(int atom_idx) const;

    /// q-point source: a q list file overrides the Monkhorst-Pack q mesh
    void set_qfile(const std::string& filename);

    void set_qmesh(int nqx, int nqy, int nqz);

    void set_conv_thr(double thr);

    void set_max_iter(int max_iter);

    void set_mix_beta(double beta);

    /// acoustic sum rule on the Gamma dynamical matrix ("no" | "simple" |
    /// "crystal"); applied to the raw matrix before the LO-TO term
    void set_asr(const std::string& asr_type);
    /// response SCF mixing scheme: "rho_plain" (beta only), "rho_kerker"
    /// (screen f_g = |G+q|^2/(|G+q|^2+a^2), a^2 via set_kerker_a2),
    /// "rho_broyden" (Modified Broyden on the density), or the potential
    /// mixes "v_sc_plain" / "v_sc_broyden" (handled in solve_displacement)
    void set_mix_type(const std::string& mix_type);

    /// Kerker screening a^2 in 1/lat0^2 units (rho_kerker scheme only)
    void set_kerker_a2(double a2);

    /// Broyden history depth (broyden scheme only)
    void set_broyden_ndim(int ndim);

    /// q = 0 response switches (epsilon_inf / Born charges / LO-TO)
    void set_compute_q0(bool flag);

    void set_loto(bool flag);

    /// q->0 direction of the non-analytic (LO-TO) term: any non-null vector
    /// is normalized to a unit direction by the data layer; the default is
    /// the isotropic (1,1,1)/sqrt(3). Consumed by run() so no direction is
    /// hardcoded in the driver anymore.
    void set_loto_dir(const ModuleBase::Vector3<double>& dir);

    /// number of q points of the current list (q file or q mesh)
    int get_nq() const;

    /// direct (reciprocal-lattice fractional) coordinates of q_idx
    ModuleBase::Vector3<double> get_qvec(int q_idx) const;

    ModuleBase::Vector3<double> get_loto_dir() const;

    /// signed Gamma frequencies (cm^-1) after the LO-TO correction; empty
    /// when loto is off or the correction has not run
    std::vector<double> get_phon_freq_loto() const;

    /// formatted per-q frequency report (deterministic layout, pinned by
    /// the serial format regression test); the LO-TO variant is empty when
    /// no corrected frequencies are available
    std::string format_q_report(int q_idx) const;

    std::string format_loto_report() const;

    /// read-only access to the data layer (consumed by the dyn file
    /// writer in the esolver post-processing)
    const DFPT_PW_Data& get_data() const;

    /// the MPI rank passed at init (0 for serial); consumed by the dyn
    /// file writer for the rank-0 output gate
    int get_myrank() const;

  private:
    std::unique_ptr<Impl> pimpl_;
};

} // namespace ModuleDFPT

#endif // DFPT_PW_H
