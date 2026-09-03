#ifndef DFTU_BASE_TOOLS_H
#define DFTU_BASE_TOOLS_H

#include <complex>
#include <vector>
#include "source_base/matrix.h"

class UnitCell;
class OccupationMatrix;

/// Free functions for DFT+U PW basis calculations.
///
/// These functions are pure (no access to Plus_U_Base members) so they can be
/// unit-tested directly by including this header. The member functions in
/// dftu_base_occ.cpp call them after computing per-atom offsets and fetching
/// the relevant member state (occ_mat, pot_uterm_pw, u_current, etc.).
namespace DFTU_BASE {

/// transform pot_onsite from Pauli basis to spin basis (in-place, nspin==4 only).
///
/// pot_onsite points to the per-atom pot_onsite block of size 4 * m_size * m_size, storing
/// 4 contiguous Pauli blocks (is=0 charge, is=1 sigma_x, is=2 sigma_y,
/// is=3 sigma_z). After the transform, the same memory holds the spin
/// representation:
///   pot_onsite[0]        <- 0.5 * (pot_onsite_pauli[0] + pot_onsite_pauli[3])
///   pot_onsite[3*size]   <- 0.5 * (pot_onsite_pauli[0] - pot_onsite_pauli[3])
///   pot_onsite[size]     <- 0.5 * (pot_onsite_pauli[1] + i * pot_onsite_pauli[2])
///   pot_onsite[2*size]   <- 0.5 * (pot_onsite_pauli[1] - i * pot_onsite_pauli[2])
void pauli_to_spin_basis(std::complex<double>* pot_onsite, int m_size);

/// compute pot_onsite and energy contribution for one atom (nspin==4, spinor).
///
/// Writes 4 Pauli blocks of pot_onsite into pot_onsite (size 4 * m_size * m_size) and
/// returns the energy_u increment. Internally calls pauli_to_spin_basis
/// to convert pot_onsite to spin basis in-place.
///
/// pot_onsite:  pointer to pot_uterm_pw[pot_uterm_pw_index[iat]]
/// occ: pointer to occ_mat[iat][target_l][0][0].c (4 Pauli blocks packed)
double compute_pot_onsite_spinor(
    std::complex<double>* pot_onsite,
    const double* occ,
    double u_value,
    double diag_coeff,
    double weight_eu,
    int m_size);

/// compute pot_onsite and energy contribution for one atom, one spin channel
/// (nspin==1 or nspin==2). Returns the energy_u increment.
///
/// pot_onsite:  pointer to the spin channel's pot_onsite block (size m_size * m_size)
/// occ: pointer to occ_mat[iat][target_l][0][is].c for this channel
double compute_pot_onsite_scalar(
    std::complex<double>* pot_onsite,
    const double* occ,
    double u_value,
    double diag_coeff,
    double weight_eu,
    int m_size);

/// accumulate occ_mat from becp for one atom, one k-point (nspin==4, spinor).
///
/// occ_mat_out points to occ_mat[iat][target_l][0][0].c, which packs 4
/// Pauli blocks contiguously (each of size tlp1*tlp1). The function adds
/// the contributions from all nbands bands for the given k-point.
///
/// becp:     projector-bra overlap for this k-point
/// npol:     2 for spinor (nspin==4)
/// nkb:      number of projectors per band per spin
/// begin_ih: offset of this atom's projectors in becp
/// m_begin:  offset of the correlated l's first m within the atom's projectors
/// tlp1:     2*target_l + 1
void accumulate_occ_spinor(
    double* occ_mat_out,
    const std::complex<double>* becp,
    int nbands,
    int npol,
    int nkb,
    int begin_ih,
    int m_begin,
    int tlp1,
    const ModuleBase::matrix& wg,
    int ik);

/// accumulate occ_mat from becp for one atom, one k-point (nspin==1 or 2).
///
/// occ_mat_out points to occ_mat[iat][target_l][0][is].c, a single channel
/// of size tlp1*tlp1. The caller selects the spin channel by passing the
/// corresponding occ_mat pointer; this function does not need is.
/// Adds contributions from all nbands bands.
///
/// becp, nbands, nkb, begin_ih, m_begin, tlp1: same as accumulate_occ_spinor
void accumulate_occ_scalar(
    double* occ_mat_out,
    const std::complex<double>* becp,
    int nbands,
    int nkb,
    int begin_ih,
    int m_begin,
    int tlp1,
    const ModuleBase::matrix& wg,
    int ik);

/// reduce occ_mat across all k-pools (per-atom, nspin-aware).
///
/// Each k-pool only accumulates occ_mat contributions from the k-points it
/// owns; this sums them across pools so occmat holds the full result.
/// nspin=1: single channel, size elements
/// nspin=2: two channels (spin-up/down) reduced separately
/// nspin=4: 4 Pauli blocks packed contiguously, reduced in one shot
void reduce_occ_mat(const UnitCell& cell,
                    const int nspin,
                    const int kpar,
                    const std::vector<int>& orbital_corr,
                    OccupationMatrix& occmat);

/// compute effective potential pot_onsite and DFT+U energy from occ_mat.
///
/// Preconditions:
///   - occmat has been accumulated from psi and reduced across k-pools.
///
/// Outputs:
///   - pot_uterm_pw: pot_onsite = U * (diag*delta - occ) written per atom
///     nspin=4: 4 Pauli blocks per atom, then transformed to spin basis
///     nspin=1: single channel
///     nspin=2: two channels in split layout [all_up | all_dn]
///   - energy_u (out): E_U = sum U * weight_eu * occ(m2,m1) * occ(m1,m2),
///     overwritten with the total energy of this call
void compute_pot_uterm_and_energy(const UnitCell& cell,
                                  const int nspin,
                                  const std::vector<double>& u_current,
                                  const std::vector<int>& orbital_corr,
                                  const std::vector<int>& pot_uterm_pw_index,
                                  const OccupationMatrix& occmat,
                                  std::vector<std::complex<double>>& pot_uterm_pw,
                                  double& energy_u);

/// accumulate occ_mat from psi for all k-points (per-device template).
///
/// Explicitly instantiated for DEVICE_CPU (and DEVICE_GPU when available)
/// in dftu_base_occ.cpp.
template <typename Device>
void accumulate_occ_one_k(const void* psi_in,
                          const ModuleBase::matrix& wg_in,
                          const UnitCell& cell,
                          const int* isk,
                          const int nspin,
                          const std::vector<int>& orbital_corr,
                          OccupationMatrix& occmat);

} // namespace DFTU_BASE

#endif
