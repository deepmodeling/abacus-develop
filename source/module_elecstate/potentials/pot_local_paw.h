#ifndef POTLOCALPAW_H
#define POTLOCALPAW_H

#include "module_base/matrix.h"
#include "pot_base.h"

namespace elecstate
{

class PotLocal_PAW : public PotBase
{
  public:
    PotLocal_PAW(){};

    void cal_fixed_v(double* vl_pseudo) override;
};

} // namespace elecstate

#endif