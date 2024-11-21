#include "lr_util.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_parameter/parameter.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"

// Here will be  the only place where GlobalCs are used (to be moved) in module_lr
#include "module_hamilt_pw/hamilt_pwdft/global.h"

template <typename T>
using D2 = void(*) (T**, size_t);
// template <typename T>
// using D3 = void(*) (T***, size_t, size_t);
// template <typename T>
// D2<T> d2 = LR_Util::_deallocate_2order_nested_ptr<T>;
// template <typename T>
// D3<T> d3 = LR_Util::delete_p3<T>;
// Change to C++ 11
D2<double> d2 = LR_Util::_deallocate_2order_nested_ptr<double>;
// D3<double> d3 = LR_Util::delete_p3<double>;


Gint& Gint::operator=(Gint&& rhs)
{
    if (this == &rhs) {return *this;
}

    this->nbx = rhs.nbx;
    this->nby = rhs.nby;
    this->nbz = rhs.nbz;
    this->ncxyz = rhs.ncxyz;
    this->nbz_start = rhs.nbz_start;
    this->bx = rhs.bx;
    this->by = rhs.by;
    this->bz = rhs.bz;
    this->bxyz = rhs.bxyz;
    this->nbxx = rhs.nbxx;
    this->ny = rhs.ny;
    this->nplane = rhs.nplane;
    this->startz_current = rhs.startz_current;

    this->gridt = rhs.gridt;
    this->ucell = rhs.ucell;

    // move hR after refactor
    this->hRGint = rhs.hRGint;
    rhs.hRGint = nullptr;
    this->hRGintCd = rhs.hRGintCd;
    rhs.hRGintCd = nullptr;
    for (int i = 0; i < this->DMRGint.size(); i++)
    {
        delete this->DMRGint[i];
    }
    for (int i = 0; i < this->hRGint_tmp.size(); i++)
    {
        delete this->hRGint_tmp[i];
    }
    this->pvdpRx_reduced = std::move(rhs.pvdpRx_reduced);
    this->pvdpRy_reduced = std::move(rhs.pvdpRy_reduced);
    this->pvdpRz_reduced = std::move(rhs.pvdpRz_reduced);
    this->DMRGint = std::move(rhs.DMRGint);
    this->hRGint_tmp = std::move(rhs.hRGint_tmp);
    this->DMRGint_full = rhs.DMRGint_full;
    rhs.DMRGint_full = nullptr;

    return *this;
}

Gint_Gamma& Gint_Gamma::operator=(Gint_Gamma&& rhs)
{
    if (this == &rhs) {return *this;
}
    Gint::operator=(std::move(rhs));

    // DM may not needed in beyond DFT ESolver
    // if (this->DM != nullptr) d3<double>(this->DM, PARAM.inp.nspin, gridt.lgd);
    assert(this->DM == nullptr);
    return *this;
}

Gint_k& Gint_k::operator=(Gint_k&& rhs)
{
    if (this == &rhs) {return *this;
}
    this->Gint::operator=(std::move(rhs));
    return *this;
}