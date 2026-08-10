#pragma once
#include "source_base/tool_title.h"
#include "source_base/ylm.h"
#include "source_basis/module_ao/orb_read.h"
#include "source_cell/atom_spec.h"
#include "source_cell/klist.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_io/module_output/cube_io.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_lcao/module_lr/dm_trans/dm_trans.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"

#include <array>
#include <complex>
#include <string>
#include <vector>

// ============================================================================
// OrbitalEvaluator: single-point NAO wavefunction evaluation
// Replicates the cubic Hermite interpolation from GintAtom::set_phi()
// ============================================================================
struct OrbitalEvaluator
{
    struct AtomTypeData
    {
        struct RadialBlock
        {
            int begin_iw = 0;
            int size = 0;
            int ylm_begin = 0;
            const double* psi_uniform = nullptr;
            const double* dpsi_uniform = nullptr;
        };

        int nw = 0, nwl = 0;                 ///< number of NAOs, max L for this atom type
        double rcut = 0.0, dr_uniform = 0.0; ///< orbital cutoff radius, uniform grid spacing
        std::vector<RadialBlock> radial_blocks;
    };
    std::vector<AtomTypeData> atypes_;

    /// @brief Initialize evaluator by extracting radial function pointers from LCAO_Orbitals
    /// @param orb numerical atomic orbital data (one Numerical_Orbital per atom type)
    void init(const LCAO_Orbitals& orb, const UnitCell& ucell);

    /// @brief Evaluate all NAO values for one atom type at relative position dr (in Bohr)
    /// @param it atom type index
    /// @param dr Cartesian vector from atom center to evaluation point (Bohr)
    /// @param phi_out array of length nw[it], filled with phi values
    void eval_phi(int it, const ModuleBase::Vector3<double>& dr, double* phi_out) const;

    /// @brief Evaluate Bloch-summed NAO values including neighboring cell contributions
    /// phi^B_mu(r, k) = Sum_R exp(i * 2pi * k_d·R_int) * phi_mu(r - tau - R).
    /// The R-sum runs over all lattice vectors within the orbital cutoff radius.
    /// Output is complex because the Bloch phase contributes a complex factor.
    /// The caller must zero-initialize phi_bloch before calling (this function accumulates).
    /// @param r evaluation point in Bohr
    /// @param ucell unit cell
    /// @param phi_bloch output array of length total_naos (complex), values are added
    /// @param kvec_d direct-coordinate k-vector
    void eval_phi_all_bloch(const ModuleBase::Vector3<double>& r,
                            const UnitCell& ucell,
                            std::complex<double>* phi_bloch,
                            const ModuleBase::Vector3<double>& kvec_d) const;

    /// @brief Evaluate a Bloch-summed LCAO KS wavefunction at r_fix point.
    /// Uses the same neighboring-cell AO sum as eval_phi_all_bloch().
    template <typename T>
    std::complex<double> eval_wfc_bloch(const ModuleBase::Vector3<double>& r_fix,
                                        int ik,
                                        int iband,
                                        const psi::Psi<T>& psi_ks,
                                        const UnitCell& ucell,
                                        const ModuleBase::Vector3<double>& kvec_d) const
    {
        psi_ks.fix_k(ik);
        int total_naos = 0;
        for (int iat = 0; iat < ucell.nat; ++iat)
        {
            total_naos += atypes_[ucell.iat2it[iat]].nw;
        }

        std::vector<std::complex<double>> phi_bloch(total_naos, std::complex<double>(0.0, 0.0));
        eval_phi_all_bloch(r_fix, ucell, phi_bloch.data(), kvec_d);

        std::complex<double> result(0.0, 0.0);
        for (int iw = 0; iw < total_naos; ++iw)
        {
            result += std::complex<double>(psi_ks(iband, iw)) * phi_bloch[iw];
        }
        return result;
    }
};

namespace LR_Util
{
template <typename T>
class ExcitonPlotter
{
    /// @brief Construct an exciton density plotter for TDA excitons.
    /// Supports average density (integrating out one particle coordinate) and
    /// conditional density (fixing the hole position and evaluating the coherent
    /// electron wavefunction).
    /// @param nspin_global global number of spin channels
    /// @param naos number of atomic-orbital basis functions
    /// @param nocc number of occupied bands per spin
    /// @param nvirt number of virtual bands per spin
    /// @param psi_ks_in Kohn-Sham wavefunction coefficients
    /// @param ucell_in unit cell
    /// @param kv_in k-point list
    /// @param gd_in grid driver for neighbor-list search
    /// @param orb_cutoff_in orbital cutoff radii per atom type
    /// @param Pgrid_in parallel grid for real-space decomposition
    /// @param rho_basis_in plane-wave basis for the real-space grid
    /// @param pX_in 2D block-cyclic distribution for X matrix
    /// @param pc_in 2D block-cyclic distribution for KS orbitals (reserved for distributed plotting)
    /// @param pmat_in parallel distribution for AO matrices
    /// @param output_dir_in output directory for generated density files
    /// @param two_fermi_in whether spin channels use separate Fermi energies
    /// @param eig excitation energies in Ry
    /// @param X BSE eigenvector amplitudes (TDA)
    /// @param openshell true for spin-unrestricted (nspin_x=2)
    /// @param orb optional LCAO_Orbitals pointer (required for conditional density)
  public:
    ExcitonPlotter(const int nspin_global,
                   const int naos,
                   const std::vector<int>& nocc,
                   const std::vector<int>& nvirt,
                   const psi::Psi<T>& psi_ks_in,
                   const UnitCell& ucell_in,
                   const K_Vectors& kv_in,
                   const Grid_Driver& gd_in,
                   const std::vector<double>& orb_cutoff_in,
                   const Parallel_Grid& Pgrid_in,
                   const ModulePW::PW_Basis& rho_basis_in,
                   const std::vector<Parallel_2D>& pX_in,
                   const Parallel_2D& pc_in,
                   const Parallel_Orbitals& pmat_in,
                   const std::string& output_dir_in,
                   const double* eig,
                   const T* X,
                   const bool openshell,
                   const LCAO_Orbitals* orb)
        : nspin_x(openshell ? 2 : 1), naos(naos), nocc(nocc), nvirt(nvirt),
          nk(nspin_global == 2 ? kv_in.get_nks() / 2 : kv_in.get_nks()),
          ldim(nk * (nspin_x == 2 ? pX_in[0].get_local_size() + pX_in[1].get_local_size() : pX_in[0].get_local_size())),
          eig(eig), X(X), kv(kv_in), ucell(ucell_in), orb_cutoff_(orb_cutoff_in), gd_(gd_in), Pgrid(Pgrid_in),
          rho_basis(rho_basis_in), output_dir_(output_dir_in), pc(pc_in), pmat(pmat_in), orb_(orb)
    {
        for (int is = 0; is < nspin_global; ++is)
        {
            psi_ks_vec.emplace_back(LR_Util::get_psi_spin(psi_ks_in, is, nk));
        }
        if (orb)
        {
            orb_eval_.init(*orb, ucell_in);
        }
    };
    // ---------------------------------------------------------------------
    // The average density integrates out one particle coordinate.
    // The conditional density fixes the hole position and evaluates the coherent electron wavefunction.
    //
    // Both use Hermitian effective DMK → DMR → cal_gint_rho (average) or
    // direct coherent real-space summation (conditional) to produce proper
    // non-negative densities.

    /// @brief Compute effective DMK for average hole density
    /// Produces D_hole(k)^T, where D_hole(k) = C_occ(k) * (X_k^H * X_k) * C_occ(k)^H,
    /// matching the AO-index order expected by DensityMatrix::cal_DMR().
    /// Marginalizes over conduction bands: M_k[v,v'] = Sum_c A_{kvc} * conj(A_{kv'c}).
    /// @param istate BSE state index
    /// @return dmk_per_kpoint
    std::vector<container::Tensor> cal_effective_dmk_hole(const int istate);

    /// @brief Compute effective DMK for average electron density
    /// Produces D_elec(k)^T, where D_elec(k) = C_virt(k) * (X_k * X_k^H) * C_virt(k)^H,
    /// matching the AO-index order expected by DensityMatrix::cal_DMR().
    /// Marginalizes over valence bands: N_k[c,c'] = Sum_v A_{kvc} * conj(A_{kvc'}).
    /// @param istate BSE state index
    /// @return dmk_per_kpoint
    std::vector<container::Tensor> cal_effective_dmk_elec(const int istate);

    /// @brief Plot average hole or electron density as a .cube file
    /// Integrates out the other particle coordinate. The effective DMK is Hermitian,
    /// so DMR's real part suffices for cal_gint_rho — no squaring needed.
    /// @param istate BSE state index
    /// @param type "hole" for average hole density, "elec" for average electron density
    void plot_average_density(const int istate, const std::string& type);

    /// @brief Plot average hole or electron density on a 2D cross-section
    /// Evaluates the Hermitian effective density matrix directly on the requested
    /// plane and writes a data file readable by plot_cond_silce.py.
    /// @param istate BSE state index
    /// @param type "hole" for average hole density, "elec" for average electron density
    /// @param plane cross-section plane: "ab", "bc", or "ca"
    /// @param slice_pos offset along the perpendicular direction (Bohr)
    /// @param npoints desired grid resolution
    /// @param range in-plane cell range {ustart, uend, vstart, vend}
    void plot_average_slice(const int istate,
                            const std::string& type,
                            const std::string& plane,
                            double slice_pos,
                            int npoints,
                            const std::vector<int>& range);

    /// @brief Compute conditional electron or hole density on a 2D BvK cross-section
    /// Evaluates psi_cond(r) on a regular grid in the chosen plane (ab, bc, or ca),
    /// spanning an explicitly selected range of primitive cells. Uses full
    /// Bloch-summed orbital evaluation with home-cell caching for efficiency.
    /// Writes a data file with grid + metadata readable by plot_cond_silce.py.
    /// @param istate BSE state index
    /// @param r_fix fixed opposite-particle position {x, y, z} in Bohr
    /// @param plane cross-section plane: "ab", "bc", or "ca"
    /// @param slice_pos offset along the perpendicular direction (Bohr)
    /// @param npoints desired grid resolution (pts per cell ≈ npoints / total_cells)
    /// @param range in-plane cell range {ustart, uend, vstart, vend}
    /// @param type "elec" or "hole"
    void plot_cond_slice(const int istate,
                         const std::array<double, 3>& r_fix,
                         const std::string& plane,
                         double slice_pos,
                         int npoints,
                         const std::vector<int>& range,
                         const std::string& type);

  private:
    /// @brief Build one effective AO coefficient vector per k-point for a
    /// conditional electron or hole wavefunction.
    std::vector<std::vector<std::complex<double>>> build_conditional_coefficients(int istate,
                                                                                  const std::array<double, 3>& r_fix,
                                                                                  bool plot_electron);

    const int nspin_x = 1;                  ///< 1 for singlet/triplet, 2 for up/down (open-shell)
    const int naos;                         ///< number of atomic-orbital basis functions
    const std::vector<int>& nocc;           ///< number of occupied bands per spin
    const std::vector<int>& nvirt;          ///< number of virtual bands per spin
    const int nk;                           ///< number of k-points
    const int ldim;                         ///< local leading dimension of X (data size per state)
    const double* eig;                      ///< excitation energies (Ry)
    const T* X;                             ///< BSE eigenvector amplitudes
    const K_Vectors& kv;                    ///< k-point list
    std::vector<psi::Psi<T>> psi_ks_vec;    ///< KS wavefunction per spin
    const UnitCell& ucell;                  ///< unit cell
    const std::vector<double>& orb_cutoff_; ///< orbital cutoff radii
    const Grid_Driver& gd_;                 ///< grid driver for neighbor search
    const Parallel_Grid& Pgrid;             ///< parallel real-space grid
    const ModulePW::PW_Basis& rho_basis;    ///< plane-wave basis for density grid
    const std::string output_dir_;          ///< output directory for density files

    const Parallel_2D& pc;               ///< KS-orbital distribution reserved for MPI plotting
    const Parallel_Orbitals& pmat;       ///< parallel distribution for AO matrices
    const LCAO_Orbitals* orb_ = nullptr; ///< orbital data for single-point evaluation
    OrbitalEvaluator orb_eval_;          ///< single-point NAO evaluator (initialized if orb_ != nullptr)
};

} // namespace LR_Util
