#ifndef DFTU_TOOLS_PW_H
#define DFTU_TOOLS_PW_H

#include <complex>
#include "source_base/matrix.h"

/// Free functions for DFT+U PW basis calculations.
///
/// These functions are pure (no access to Plus_U_Base members) so they can be
/// unit-tested directly by including this header. The member functions in
/// dftu_pw.cpp call them after computing per-atom offsets and fetching the
/// relevant member state (occ_mat, pot_uterm_pw, u_current, etc.).
namespace dftu_pw {

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

} // namespace dftu_pw

#endif
