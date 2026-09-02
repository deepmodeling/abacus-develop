#ifndef DFTU_YUKAWA_H
#define DFTU_YUKAWA_H

class Plus_U;
class UnitCell;

#ifdef __LCAO
namespace DFTU_LCAO {

/**
 * @brief Spherical modified Bessel function of the first kind (Yukawa kernel).
 *
 * Evaluates i_k(r*lambda) for even orders k = 0, 2, 4, 6 used in the Slater
 * integral construction for the Yukawa-screened DFT+U. The implementation uses
 * piecewise small-x expansions and large-x closed forms, so the caller must
 * pass an order handled by those branches; other orders return 0.
 *
 * @param k      order (supported: 0, 2, 4, 6)
 * @param r      radial distance
 * @param lambda Yukawa screening length
 * @return       function value
 */
double spherical_Bessel(const int k, const double r, const double lambda);

/**
 * @brief Spherical modified Hankel function of the second kind (Yukawa kernel).
 *
 * Evaluates k_k(r*lambda) for even orders k = 0, 2, 4, 6 used together with
 * spherical_Bessel in the radial Slater integrals. Same piecewise scheme as
 * spherical_Bessel; unsupported orders return 0.
 *
 * @param k      order (supported: 0, 2, 4, 6)
 * @param r      radial distance
 * @param lambda Yukawa screening length
 * @return       function value
 */
double spherical_Hankel(const int k, const double r, const double lambda);

/**
 * @brief Determine the Yukawa screening length lambda and store it on dftu.
 *
 * If dftu carries a positive yukawa_lambda configuration value, that value is
 * used directly; otherwise lambda is estimated from the charge density
 * (Thomas-Fermi-like) and rescaled by 1.6. The spin channel count is read
 * from the global PARAM.inp.nspin rather than a Plus_U member.
 *
 * @param dftu Plus_U state (lambda is written back via set_lambda)
 * @param rho  charge density per spin channel
 * @param nrxx number of real-space grid points
 */
void cal_yukawa_lambda(Plus_U& dftu, double** rho, const int& nrxx);

/**
 * @brief Compute the Slater integrals Fk for the correlated orbital of atom type T.
 *
 * Accumulates into dftu.Fk[T][L][chi][k] using the radial orbital grid and the
 * Yukawa-screened spherical Bessel/Hankel kernels. Only acts when Yukawa
 * screening is enabled on dftu.
 *
 * @param dftu  Plus_U state (provides ptr_orb, lambda, Fk, use_yukawa)
 * @param ucell unit cell
 * @param L     angular momentum
 * @param T     atom type
 */
void cal_slater_Fk(Plus_U& dftu, const UnitCell& ucell, const int L, const int T);

/**
 * @brief Compute Yukawa-screened Slater integrals and derive U/J for all atoms.
 *
 * Drives cal_yukawa_lambda then cal_slater_Fk over correlated orbitals and
 * writes the resulting U_Yukawa/J_Yukawa/u_current back onto dftu. No-op when
 * Yukawa screening is disabled on dftu.
 *
 * @param dftu Plus_U state (U_Yukawa, J_Yukawa, u_current, lambda written back)
 * @param ucell unit cell
 * @param rho   charge density per spin channel
 * @param nrxx  number of real-space grid points
 */
void cal_slater_UJ(Plus_U& dftu, const UnitCell& ucell, double** rho, const int& nrxx);

} // namespace DFTU_LCAO
#endif

#endif
