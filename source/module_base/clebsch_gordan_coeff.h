#ifndef CLEBSCH_GORDAN_H
#define CLEBSCH_GORDAN_H

#include "module_base/intarray.h"
#include "module_base/matrix.h"
#include "module_base/realarray.h"
#include "module_base/vector3.h"

namespace ModuleBase
{

class Clebsch_Gordan
{
  public:
    Clebsch_Gordan();
    ~Clebsch_Gordan();

    /**
     * @brief computes Clebsch-Gordan coefficient
     *
     * this routine computes the coefficients of the expansion of the product
     * of two real spherical harmonics into real spherical harmonics.
     *
     *    Y_limi(r) * Y_ljmj(r) = \sum_LM  ap(LM,limi,ljmj)  Y_LM(r)
     *
     * The indices limi,ljmj and LM assume the order for real spherical
     * harmonics given in routine ylmr2
     *
     * @param lli [in] lmaxkb + 1: the maximum li considered
     * @param ap [out] coefficients of the expansion
     * @param lpx [out] for each input limi,ljmj is the number of LM in the sum
     * @param lpl [out] for each input limi,ljmj points to the allowed LM
     */
    static void clebsch_gordan(const int& lli,
                               ModuleBase::realArray& ap,
                               ModuleBase::IntArray& lpx,
                               ModuleBase::IntArray& lpl);

  private:
    /**
     * @brief generate random vector
     *
     * @param llx [in] the number of vectors
     * @param r [out] an array of vectors
     */
    static void gen_rndm_r(const int& llx, ModuleBase::Vector3<double>* r);

    /**
     * @brief store the inverse of ylm(ir,lm) in mly(lm,ir)
     *
     * @param L [in] angular momentum L
     * @param li [in] angular momentum li
     * @param lj [in] angular momentum lj
     * @param llx [in] the number of vectors
     * @param ylm [in] real spherical harmonics
     * @param mly [in] the inverse of ylm(ir,lm)
     * @return double the expansion coefficients
     */
    static double compute_ap(const int& L,
                             const int& li,
                             const int& lj,
                             const int& llx,
                             const ModuleBase::matrix& ylm,
                             const ModuleBase::matrix& mly);
};

} // namespace ModuleBase

#endif // CLEBSCH_GORDAN_H