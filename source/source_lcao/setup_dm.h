#ifndef SETUP_DM_H
#define SETUP_DM_H 

#include "source_estate/module_dm/density_matrix.h"

#include <vector>

namespace LCAO_DOMAIN
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

    // allocate density matrix
    void allocate_dm(const K_Vectors* kv, const Parallel_Orbitals* paraV, const int nspin);

    elecstate::DensityMatrix<TK, double>* get_dm() const
    {
        return const_cast<elecstate::DensityMatrix<TK, double>*>(this->dm);
    }

    private:

    elecstate::DensityMatrix<TK, double>* dm = nullptr;

};


} // namespace elecstate

#endif
