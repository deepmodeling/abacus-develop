#ifndef MODULE_DM_DM_HOLDER_H
#define MODULE_DM_DM_HOLDER_H

#include "source_estate/module_dm/density_matrix.h"

namespace module_dm
{
template <typename TK>
class Setup_DM
{
  public:
    Setup_DM()
    {
    } // will be called by ElecStateLCAO_TDDFT

    ~Setup_DM()
    {
        if (this->dm != nullptr)
        {
            delete this->dm;
        }
    }

    DensityMatrix<TK, double>* dm = nullptr;
};

} // namespace module_dm

#endif
