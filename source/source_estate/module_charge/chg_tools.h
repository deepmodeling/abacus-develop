#ifndef CHG_TOOLS_H
#define CHG_TOOLS_H

// Free numerical kernels extracted from the Charge class so that the
// density math (summation, electron-count integration, non-linear core
// correction) can be tested and reused without dragging in Charge's state.
// All grid / geometry inputs are passed explicitly instead of being read
// from Charge members or globals.

#include <complex>
#include <vector>

class UnitCell;
namespace ModuleBase
{
class ComplexMatrix;
template <class T>
class Vector3;
}
namespace ModulePW
{
class PW_Basis;
}

namespace module_charge
{

// Compute the core charge (non-linear core correction) on the real-space
// 3D mesh. rho_core / rhog_core are the output buffers previously owned by
// Charge; rhopw supplies the reciprocal-grid geometry and FFT backend.
void set_rho_core(const UnitCell& ucell,
                  const ModuleBase::ComplexMatrix& structure_factor,
                  const bool* numeric,
                  double* rho_core,
                  std::complex<double>* rhog_core,
                  const ModulePW::PW_Basis& rhopw);

// First-order core-density response (the q-shifted derivative of the
// core charge, QE addcore convention). Conventions: q_cart is in the
// same scaled Cartesian units as rhopw.gcar (q_cart = q_frac * ucell.G),
// so the physical wavevector is K = tpiba * (gcar + q_cart) in 1/bohr;
// with tau_cart = lat0 * tau (bohr), for the displacement of atom iat
// along direction dir,
//   drhoc_g(G) = -i K_dir * rho_atc~(|K|) * exp(-i K.tau_cart)
// where rho_atc~ is the numeric Simpson radial transform
//   4pi/Omega Int r^2 rho_atc(r) j0(|K| r) dr.
// Unlike set_rho_core this evaluates per G (|K| differs at every ig, so
// the gg_uniq shell cache cannot be reused) and carries the
// single-atom structure-factor phase. Empty output for a non-nlcc atom
// or an out-of-range (iat, dir); the |K| = 0 harmonic is left zero.
// Consumed by the DFPT first-order density (DFPT_Rho::compute_drhoc_g
// delegates here).
void set_drho_core(const UnitCell& ucell,
                   int iat,
                   int dir,
                   const ModuleBase::Vector3<double>& q_cart,
                   const ModulePW::PW_Basis& rhopw,
                   std::vector<std::complex<double>>& drhoc_g);

// Core-charge form factor at a single physical wavevector k (1/bohr):
//   4pi/Omega * Int r^2 rho_atc(r) j0(k r) dr   (numeric Simpson).
// The shared radial kernel of every per-point core-charge evaluation
// (the first-order response above and the DFPT dynmatcc second-order
// route); the ground-state set_rho_core instead shells it through the
// cached, multithreaded non_linear_core_correction. scratch is a
// caller-owned buffer of size mesh reused across calls (the kernel
// allocates nothing per G).
double core_charge_form_factor(int mesh,
                               const double* r,
                               const double* rab,
                               const double* rho_atc,
                               double omega,
                               double k,
                               std::vector<double>& scratch);

// Sum the (spin-resolved) charge density over the real-space grid and
// convert it to a total charge using the cell volume and grid size.
// rho[is][ir] is the density for spin channel is. nspin0 is the number of
// spin channels to include (2 for nspin==2, 1 otherwise).
double sum_rho(double* const* rho,
               const int nspin0,
               const int nrxx,
               const double omega,
               const int nxyz);

// Integrate a single spin channel rho_in over the grid to obtain the
// electron number, scaled by omega / nxyz. Reduction over the pool is
// performed internally under __MPI.
double cal_rho2ne(const double* rho_in,
                  const int nrxx,
                  const double omega,
                  const int nxyz);

// Check that the charge density integrates to the target electron number.
// rho[is][ir] is the density for spin channel is. For nspin == 1 or 4 the
// total density in rho[0] is compared with nelec; for nspin == 2 the spin-up
// (rho[0]) and spin-down (rho[1]) channels must each integrate to a
// non-negative number and their sum must match nelec. A mismatch emits a
// warning; a negative spin-channel electron number aborts the run.
void check_rho(double* const* rho,
               const int nspin,
               const int nrxx,
               const double omega,
               const int nxyz,
               const double nelec);

// Non-linear core correction: Fourier transform of the (numeric) core
// charge. gg_uniq / ngg supply the reciprocal grid shells previously read
// from Charge::rhopw.
/// Radial grid and density data for non-linear core correction
struct NlcCtx
{
    bool numeric;       ///< whether the pseudo-potential has numeric data
    double omega;      ///< cell volume
    double tpiba2;     ///< 2*pi/alat squared
    int mesh;          ///< radial mesh size
    const double* r;   ///< radial grid points
    const double* rab; ///< radial grid weights
    const double* rhoc; ///< core charge on the radial grid
    const double* gg_uniq; ///< unique |G|^2 values
    int ngg;           ///< number of unique |G|^2 shells
};

void non_linear_core_correction(const NlcCtx& ctx,
                                double* rhocg);

} // namespace module_charge

#endif // CHG_TOOLS_H
