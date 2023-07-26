#ifndef TWO_CENTER_INTEGRATOR_H_
#define TWO_CENTER_INTEGRATOR_H_

#include "module_basis/module_nao/two_center_table.h"
#include "module_basis/module_nao/real_gaunt_table.h"
#include "module_basis/module_nao/radial_collection.h"
#include "module_base/vector3.h"

/**
 * @brief A class to compute two-center integrals
 *
 * This class computes two-center integrals
 *
 *                     /    
 *              I(R) = | dr phi1(r) (op) phi2(r - R)
 *                     /               
 *
 * as well as their gradients, where op is 1 (overlap) or minus Laplacian (kinetic),
 * and phi1, phi2 are "atomic-orbital-like" functions of the form
 *
 *              phi(r) = chi(|r|) * Ylm(r/|r|)
 *
 * where chi is some numerical radial function and Ylm is some real spherical harmonics.
 *
 * This class is designed to efficiently compute the two-center integrals between 
 * two "collections" of the above functions with various R, e.g., the overlap integrals 
 * between all numerical atomic orbitals and all Kleinman-Bylander nonlocal projectors,
 * the overlap & kinetic integrals between all numerical atomic orbitals, etc.
 * This is done by tabulating the radial part of the integrals on an r-space grid and
 * the real Gaunt coefficients in advance.
 *
 * See the developer's document for more details.
 *                                                                                      */
class TwoCenterIntegrator
{
  public:
    TwoCenterIntegrator();
    TwoCenterIntegrator(const TwoCenterIntegrator&) = delete;
    TwoCenterIntegrator& operator=(const TwoCenterIntegrator&) = delete;

    ~TwoCenterIntegrator();

    /**
     * @brief Tabulate the radial part of the two-center integrals
     *
     * @param[in] bra          The radial functions of the first collection.
     * @param[in] ket          The radial functions of the second collection.
     * @param[in] op           Operator, could be 'S' or 'T'.
     * @param[in] nr           Number of r-space grid points.
     * @param[in] cutoff       r-space cutoff radius.
     * @param[in] with_deriv   If true, the derivative of radial table is also tabulated.
     *                         This is necessary to compute the gradient of integrals.
     * @param[in] rgt          Pointer to a real Gaunt table.
     *                                                                                  */
    void tabulate(const RadialCollection& bra,
                  const RadialCollection& ket,
                  const char op,
                  const int nr,
                  const double cutoff,
                  const bool with_deriv,
                  RealGauntTable* const rgt = nullptr
    );

    /**
     * @brief Compute the two-center integrals
     *
     * This function calculates the two-center integral
     *
     *                     /    
     *              I(R) = | dr phi1(r) (op_) phi2(r - R)
     *                     /               
     *
     * or its gradient by using the tabulated radial part and real Gaunt coefficients.
     *
     * @param[in] itype1       Element index of orbital 1.
     * @param[in] l1           Angular momentum of orbital 1.
     * @param[in] izeta1       Zeta number of orbital 1.
     * @param[in] m1           Magnetic quantum number of orbital 1.
     * @param[in] itype2       Element index of orbital 2.
     * @param[in] l2           Angular momentum of orbital 2.
     * @param[in] izeta2       Zeta number of orbital 2.
     * @param[in] m2           Magnetic quantum number of orbital 2.
     * @param[in] vR           R2 - R1.
     * @param[in] deriv        If true, the gradient of integrals is computed.
     * @param[out] out         The computed integrals. If deriv is true, out[0], out[1]
     *                         and out[2] are the x, y, z components of the gradient.
     *                                                                                  */
    void calculate(const int itype1, 
                   const int l1, 
                   const int izeta1, 
                   const int m1, 
                   const int itype2,
                   const int l2,
                   const int izeta2,
                   const int m2,
	               const ModuleBase::Vector3<double>& vR, // vR = R2 - R1
                   const bool deriv,
                   double* out
    ) const;

  private:
    bool is_tabulated_;

    char op_;
    bool with_deriv_;
    bool use_internal_gaunt_;

    TwoCenterTable table_;
    RealGauntTable* rgt_;

    /**
     * @brief Compute the index of (l,m) in the array of spherical harmonics
     *
     * Spherical harmonics in ABACUS are stored in the following order:
     *
     * index  0   1   2   3   4   5   6   7   8   9  10  11  12 ...
     *   l    0   1   1   1   2   2   2   2   2   3   3   3   3 ...
     *   m    0   0   1  -1   0   1  -1   2  -2   0   1  -1   2 ...
     *                                                                                  */
    int ylm_index(const int l, const int m) const;
};

#endif
