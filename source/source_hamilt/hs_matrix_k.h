#ifndef HS_MATRIX_K_H
#define HS_MATRIX_K_H

#include "source_basis/module_ao/parallel_orbitals.h"

#include <complex>
#include <cstddef>
#include <vector>
namespace hamilt
{
template <typename TK>
class HS_Matrix_K
{

  public:
    HS_Matrix_K(const Parallel_Orbitals *paraV, bool no_s = false)
    {
        // size() is size_t; nrloc*ncloc can exceed INT_MAX for large
        // LCAO bases (issue #8342).
        hk.resize(paraV->para_2d().size());
        if (!no_s)
        {
            sk.resize(paraV->para_2d().size());
        }
        this->pv = paraV;
    }

    TK *get_hk() { return hk.data(); }

    TK *get_sk() { return sk.data(); }

    std::size_t get_size() const { return hk.size(); }

    void set_zero_hk() { hk.assign(hk.size(), 0); }

    void set_zero_sk() { sk.assign(sk.size(), 0); }

    const Parallel_Orbitals *get_pv() const { return this->pv; }

  private:
    std::vector<TK> hk;

    std::vector<TK> sk;

    const Parallel_Orbitals *pv = nullptr;
};
} // namespace hamilt

#endif