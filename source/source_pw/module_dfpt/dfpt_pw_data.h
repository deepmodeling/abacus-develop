// ============================================================
// This code is added by Mohan Chen on 2026-05-18.
// This code is currently in the design phase and has not been
// put into production yet. It may change in the future.
// Please use this code with caution. Only developers who know
// what they are doing should use this code.
// ============================================================

#ifndef DFPT_PW_DATA_H
#define DFPT_PW_DATA_H

#include "source_base/matrix.h"
#include "source_base/complexmatrix.h"
#include "source_base/vector3.h"
#include "source_psi/psi.h"
#include "source_cell/qlist.h"
#include <map>
#include <utility>
#include <vector>
#include <complex>

class Plus_U_Base;

namespace ModuleDFPT {

/// Occupied-band classifier shared by the projector build, the Sternheimer
/// driver, the response-density accumulation and the 2n+1 assembly. A band
/// counts as occupied iff its weight exceeds half of the per-k full
/// reference (band 0 is always the deepest, fully occupied band). A fixed
/// absolute threshold makes the Sternheimer projector jump between k
/// samplings: a smeared Fermi-tail band with weight ~1e-6 lands on either
/// side of 1e-8 depending on where the sampling's Fermi level sits, which
/// opens or closes its empty-state channel in (H-eps)^-1 and changes the
/// converged response at the percent level. The majority criterion keeps
/// that channel open for tail bands (their occupied-type contribution is
/// f-weighted and negligible), reproducing the insulator limit that
/// finite-difference references follow.
inline bool dfpt_band_occupied(const ModuleBase::matrix& wg, int ik, int ib)
{
    return wg(ik, ib) > 0.5 * wg(ik, 0);
}

class DFPT_PW_Data {
public:
    DFPT_PW_Data();
    ~DFPT_PW_Data();
    
    void init(ModuleCell::QList* qlist, int nk, int nbands, int npw_max, 
              int nrxx, int nspin, int nat, const Plus_U_Base* dftu);
    
    void clean();
    
    int get_nq() const;
    ModuleBase::Vector3<double> get_qvec(int q_idx) const;
    int get_nirr(int q_idx) const;
    std::vector<int> get_irrep_modes(int q_idx, int irrep) const;
    
    void set_dpsi(int q_idx, int k_idx, int band_idx,
                  const std::vector<std::complex<double>>& psi);
    std::vector<std::complex<double>> get_dpsi(int q_idx, int k_idx, int band_idx) const;
    
    void set_drho_r(int q_idx, int spin, const std::vector<double>& rho);
    std::vector<double> get_drho_r(int q_idx, int spin) const;
    void set_drho_g(int q_idx, int spin, const std::vector<std::complex<double>>& rho);
    std::vector<std::complex<double>> get_drho_g(int q_idx, int spin) const;
    
    void set_dv_r(int q_idx, int spin, const std::vector<double>& v);
    std::vector<double> get_dv_r(int q_idx, int spin) const;
    
    /// First-order perturbation potential dV stored as complex plane-wave
    /// coefficients (indexed by the rho-grid ig) and as the corresponding
    /// complex real-space array on the shared FFT grid (C1).
    /// The reciprocal coefficients already carry the -i(Delta+q) prefactor
    /// and the atomic phase exp(i(Delta+q).tau); the real-space array is
    /// their inverse Fourier transform, so that dV.psi is a plain cyclic
    /// convolution on the shared grid (apply_dv needs no extra q-phase).
    void set_dv_recip_c(int q_idx, int spin, const std::vector<std::complex<double>>& v);
    std::vector<std::complex<double>> get_dv_recip_c(int q_idx, int spin) const;
    void set_dv_rc(int q_idx, int spin, const std::vector<std::complex<double>>& v);
    std::vector<std::complex<double>> get_dv_rc(int q_idx, int spin) const;
    
    /// The dynamical matrix at a generic q is complex Hermitian; stored as a
    /// ModuleBase::ComplexMatrix (C5), consumed by DFPT_Phon::diagonalize
    /// through the LapackConnector::zheev wrapper.
    void set_dynmat(int q_idx, const ModuleBase::ComplexMatrix& dm);
    ModuleBase::ComplexMatrix get_dynmat(int q_idx) const;
    void set_phon_freq(int q_idx, const std::vector<double>& freq);
    std::vector<double> get_phon_freq(int q_idx) const;
    
    void set_dielectric(const ModuleBase::matrix& eps);
    ModuleBase::matrix get_dielectric() const;
    void set_born(int atom_idx, const ModuleBase::matrix& z);
    ModuleBase::matrix get_born(int atom_idx) const;
    
    void set_compute_q0(bool flag) { compute_q0_ = flag; }
    bool get_compute_q0() const { return compute_q0_; }
    void set_loto(bool flag) { loto_ = flag; }
    bool get_loto() const { return loto_; }

    /// q->0 direction of the non-analytic (LO-TO) term, as a unit vector.
    /// The setter normalizes; a null vector falls back to the isotropic
    /// default (1,1,1)/sqrt(3) (documented cubic-crystal default; a general
    /// direction control arrives with the irrep machinery of stage A).
    void set_loto_dir(const ModuleBase::Vector3<double>& dir);
    ModuleBase::Vector3<double> get_loto_dir() const { return loto_dir_; }

    /// signed Gamma frequencies (cm^-1) after the non-analytic LO-TO term
    /// along loto_dir_; empty until add_loto + diagonalize_loto have run
    void set_phon_freq_loto(const std::vector<double>& freq) { phon_freq_loto_ = freq; }
    std::vector<double> get_phon_freq_loto() const { return phon_freq_loto_; }

    /// The perturbation currently being solved: displacement of which linear
    /// atom index (over all atoms) and along which cartesian direction.
    /// Set by DFPT_Pert::build_dv and consumed by DFPT_Pert::apply_dv so the
    /// Stern solver can keep applying the same perturbation per irrep without
    /// re-passing (atom,dir) on every matrix-vector product.
    void set_pert_atom(int atom_idx) { pert_atom_ = atom_idx; }
    int get_pert_atom() const { return pert_atom_; }
    void set_pert_dir(int dir) { pert_dir_ = dir; }
    int get_pert_dir() const { return pert_dir_; }
    
    void set_is_metal(bool flag) { is_metal_ = flag; }
    bool get_is_metal() const { return is_metal_; }
    void set_dmu(double dmu) { dmu_ = dmu; }
    double get_dmu() const { return dmu_; }
    
    void set_max_iter(int iter) { max_iter_ = iter; }
    int get_max_iter() const { return max_iter_; }
    void set_conv_thr(double thr) { conv_thr_ = thr; }
    double get_conv_thr() const { return conv_thr_; }

    /// Per-(q, irrep) SCF convergence ledger (B4: sunk from the retired
    /// DFPT_IrrepData adapter). The irrep dimension is a stage-A slot:
    /// today the single fallback irrep 0 carries the full 3N displacement
    /// basis, and DFPT_PW::run records one ledger entry per outer SCF pass
    /// (worst displacement residual of the pass); missing keys read as
    /// not-converged / empty history / iteration 0.
    void set_converged(int q_idx, int irrep, bool flag);
    bool get_converged(int q_idx, int irrep) const;
    void add_residual(int q_idx, int irrep, double r);
    std::vector<double> get_residuals(int q_idx, int irrep) const;
    void set_current_iter(int q_idx, int irrep, int iter);
    int get_current_iter(int q_idx, int irrep) const;
    
    /// DFT+U interface reservation (U0):
    /// the DFPT modules never read global input state directly; the esolver
    /// layer decides whether DFT+U is active and passes a non-null
    /// Plus_U_Base* only then.
    /// with_u(): a DFT+U provider is wired (dft_plus_u enabled upstream).
    /// u_active(): the provider is additionally usable (occupation matrices
    ///             initialized, which the ground state does when DFT+U
    ///             actually runs; a provider without them stays inactive).
    bool with_u() const { return dftu_ != nullptr; }
    bool u_active() const;
    const Plus_U_Base* get_dftu() const { return dftu_; }
    
    /// first-order occupation matrix (docc) storage, indexed by q.
    /// lazy allocation: unset / out-of-range reads return an empty vector.
    void set_docc(int q_idx, const std::vector<std::complex<double>>& occ);
    std::vector<std::complex<double>> get_docc(int q_idx) const;

    /// converged screened response potential of displacement (atom, dir):
    /// the real-space q-shifted complex amplitude v_sc used by the last
    /// Sternheimer iteration of that displacement. The 2n+1 accumulation
    /// needs it to complete the term2 cross section
    /// 2<dpsi^a|dV_ext^b + dV_sc^b|psi> (screening channel).
    void set_vsc_r(int atom_idx, int dir,
                   const std::vector<std::complex<double>>& v);
    std::vector<std::complex<double>> get_vsc_r(int atom_idx, int dir) const;

    /// converged dpsi of displacement (atom, dir), indexed [k][band]; the
    /// two-pass 2n+1 accumulation reads it back after all displacements of
    /// the basis have been solved (the working dpsi slots get overwritten by
    /// later solves).
    void set_dpsi_disp(int atom_idx, int dir,
                       const std::vector<std::vector<std::vector<std::complex<double>>>>& d);
    std::vector<std::vector<std::vector<std::complex<double>>>>
    get_dpsi_disp(int atom_idx, int dir) const;

    /// conduction-projected position operator P_c r_dir |u_(k,band)> of the
    /// q = 0 mesh, solved exactly as a linear response ((H - eps_band) Y =
    /// -(i/tpiba) dH/dk_dir |u>), indexed [dir][k][band]; the screened Born
    /// charge contraction <dpsi^kappa|P_c r|u> avoids the empty-eigenvector
    /// truncation of the explicit r-matrix sum
    void set_pos_resp(int dir,
                      const std::vector<std::vector<std::vector<std::complex<double>>>>& y);
    std::vector<std::vector<std::vector<std::complex<double>>>>
    get_pos_resp(int dir) const;

    /// converged screened E-field response dpsi^E(dir) of the q = 0 mesh
    /// (QE solve_e + dfpt_kernel fixed point on the rhs
    /// -(Y^dir + dV_sc^E|psi>)), indexed [dir][k][band]
    void set_dpsi_efield(
        int dir,
        const std::vector<std::vector<std::vector<std::complex<double>>>>& d);
    std::vector<std::vector<std::vector<std::complex<double>>>>
    get_dpsi_efield(int dir) const;

private:
    ModuleCell::QList* qlist_ = nullptr;
    
    int nk_ = 0;
    int nbands_ = 0;
    int npw_max_ = 0;
    int nrxx_ = 0;
    int nspin_ = 1;
    int nat_ = 0;
    
    /// first-order wavefunction response, indexed [q][k][band]; each entry is
    /// the dpsi on the k+q basis for that band (a vector of complex coefficients).
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> dpsi_;
    
    std::vector<std::vector<std::vector<double>>> drho_r_;
    std::vector<std::vector<std::vector<std::complex<double>>>> drho_g_;
    
    std::vector<std::vector<std::vector<double>>> dv_r_;
    
    std::vector<std::vector<std::vector<std::complex<double>>>> dv_recip_c_;
    std::vector<std::vector<std::vector<std::complex<double>>>> dv_rc_;
    
    std::vector<ModuleBase::ComplexMatrix> dynmat_;
    std::vector<std::vector<double>> phon_freq_;
    
    bool compute_q0_ = false;
    bool loto_ = false;
    ModuleBase::Vector3<double> loto_dir_{1.0 / std::sqrt(3.0),
                                           1.0 / std::sqrt(3.0),
                                           1.0 / std::sqrt(3.0)};
    std::vector<double> phon_freq_loto_;
    int pert_atom_ = -1;
    int pert_dir_ = -1;
    ModuleBase::matrix dielectric_;
    std::vector<ModuleBase::matrix> born_;
    
    bool is_metal_ = false;
    double dmu_ = 0.0;
    
    /// DFT+U reservation state (U0)
    const Plus_U_Base* dftu_ = nullptr;
    std::vector<std::vector<std::complex<double>>> docc_;

    /// converged v_sc per displacement (atom, dir): [3*nat] entries
    std::vector<std::vector<std::complex<double>>> vsc_r_;

    /// converged dpsi per displacement (atom, dir): [3*nat][k][band] entries
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> dpsi_disp_;

    /// conduction-projected position response P_c r_dir|u> per direction:
    /// [3][k][band] entries
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> pos_resp_;

    /// converged E-field response dpsi^E per direction: [3][k][band]
    std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
    dpsi_efield_;
    
    int max_iter_ = 100;
    double conv_thr_ = 1e-8;

    ///< per-(q, irrep) SCF ledger (B4: absorbed from DFPT_IrrepData)
    std::map<std::pair<int, int>, bool> converged_;
    std::map<std::pair<int, int>, std::vector<double>> residuals_;
    std::map<std::pair<int, int>, int> current_iter_;
    
    bool is_initialized_ = false;
    
    void allocate_memory();
    void deallocate_memory();
};

} // namespace ModuleDFPT

#endif // DFPT_PW_DATA_H