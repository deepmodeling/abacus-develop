#include "esolver_double_xc.h"

namespace ModuleESolver
{

template <typename TK, typename TR>
ESolver_DoubleXC<TK, TR>::ESolver_DoubleXC()
{
    this->classname = "ESolver_DoubleXC";
    this->basisname = "LCAO";
}

template <typename TK, typename TR>
ESolver_DoubleXC<TK, TR>::~ESolver_DoubleXC()
{
}

template <typename TK, typename TR>
void ESolver_DoubleXC<TK, TR>::before_all_runners(UnitCell& ucell, const Input_para& inp)
{
    ModuleBase::TITLE("ESolver_DoubleXC", "before_all_runners");
    ModuleBase::timer::tick("ESolver_DoubleXC", "before_all_runners");

    ESolver_KS_LCAO<TK, TR>::before_all_runners(ucell, inp);

    ModuleBase::timer::tick("ESolver_DoubleXC", "before_all_runners");
}

template <typename TK, typename TR>
void ESolver_DoubleXC<TK, TR>::runner(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_DoubleXC", "runner");
    ModuleBase::timer::tick("ESolver_DoubleXC", "runner");

    ModuleBase::timer::tick("ESolver_DoubleXC", "runner");
}

template <typename TK, typename TR>
void ESolver_DoubleXC<TK, TR>::after_all_runners(UnitCell& ucell)
{
    ModuleBase::TITLE("ESolver_DoubleXC", "after_all_runners");
    ModuleBase::timer::tick("ESolver_DoubleXC", "after_all_runners");

    ESolver_KS_LCAO<TK, TR>::after_all_runners(ucell);

    ModuleBase::timer::tick("ESolver_DoubleXC", "after_all_runners");
};

template class ESolver_DoubleXC<double, double>;
template class ESolver_DoubleXC<std::complex<double>, double>;
template class ESolver_DoubleXC<std::complex<double>, std::complex<double>>;

} // namespace ModuleESolver
