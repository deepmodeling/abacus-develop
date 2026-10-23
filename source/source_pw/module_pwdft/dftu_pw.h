#ifndef DFTU_PW_H
#define DFTU_PW_H

#include <complex>
#include <vector>
#include "source_base/matrix.h"

class UnitCell;
class OccupationMatrix;
class OccMatMixer;
class Charge_Mixing;

/// Orchestration functions for DFT+U PW basis calculations.
///
/// These functions manage state (psi, OnsiteProjector, Charge_Mixing, etc.)
/// and call the pure functions in dftu_pw_tools.h. They cannot be unit-tested
/// directly because they depend on runtime objects.
namespace DFTU_BASE {

/// accumulate occ_mat from psi for all k-points (per-device template).
///
/// Explicitly instantiated for DEVICE_CPU (and DEVICE_GPU when available)
/// in dftu_pw.cpp.
template <typename Device>
void accumulate_occ_one_k(const void* psi_in,
                          const ModuleBase::matrix& wg_in,
                          const UnitCell& cell,
                          const int* isk,
                          const int nspin,
                          const std::vector<int>& l_channel,
                          OccupationMatrix& occmat);

/// calculate the local occupation number matrix for PW based wave functions.
///
/// This is the PW-basis entry point that:
///   1. saves and zeroes the occupation matrix
///   2. accumulates it from psi via accumulate_occ_one_k
///   3. reduces across k-pools via reduce_occ_mat
///   4. applies occupation-matrix mixing when enabled
///   5. computes the effective potential and DFT+U energy
///
/// All state is passed explicitly so this function can be unit-tested
/// without constructing a Plus_U_Base object.
void cal_occ_pw(const void* psi_in,
                const ModuleBase::matrix& wg_in,
                const UnitCell& cell,
                Charge_Mixing* p_chgmix,
                const int* isk,
                const int kpar,
                const int nspin,
                const std::string& device,
                const std::vector<int>& l_channel,
                const std::vector<double>& u_current,
                const std::vector<int>& uterm_mat_index,
                OccupationMatrix& occmat,
                OccMatMixer* occ_mixer,
                std::vector<std::complex<double>>& uterm_mat,
                double& energy_u);

} // namespace DFTU_BASE

#endif // DFTU_PW_H
