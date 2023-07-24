#ifndef TWO_CENTER_INTEGRAL_H_
#define TWO_CENTER_INTEGRAL_H_

#include "module_basis/module_nao/two_center_table.h"
#include "module_basis/module_nao/real_gaunt_table.h"
#include "module_basis/module_nao/radial_collection.h"
#include "module_base/vector3.h"

class TwoCenterIntegral
{
  public:
    TwoCenterIntegral();
    TwoCenterIntegral(const TwoCenterIntegral&) = delete;
    TwoCenterIntegral& operator=(const TwoCenterIntegral&) = delete;

    ~TwoCenterIntegral();

    void build(const RadialCollection& bra,
               const RadialCollection& ket,
               const char op,
               const int nr,
               const double cutoff,
               const bool with_deriv,
               RealGauntTable* const rgt
    );

    void get(const int itype1, 
             const int l1, 
             const int izeta1, 
             const int m1, 
             const int itype2,
             const int l2,
             const int izeta2,
             const int m2,
	         const ModuleBase::Vector3<double>& dR, // dR = R2 - R1
             const bool deriv,
             double* out
    ) const;

  private:
    bool with_deriv_;
    bool use_internal_gaunt_;

    TwoCenterTable table_;
    RealGauntTable* rgt_;

    int lm_index(const int l, const int m) const;
};

#endif
