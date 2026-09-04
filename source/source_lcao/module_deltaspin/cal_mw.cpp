#ifdef __LCAO

/**
 * @file cal_mw.cpp
 * @brief Thin LCAO shells on SpinConstrain: cal_mi_lcao() and set_operator().
 *
 * The actual LCAO magnetic-moment implementation lives in
 * deltaspin_lcao_mi.cpp as free functions over ScState; the member
 * functions below only adapt the singleton's stored pointers.
 */

#include "spin_constrain.h"

#include "deltaspin_lcao_mi.h"
#include "source_lcao/module_operator_lcao/dspin_lcao.h"
#include "source_estate/module_dm/density_matrix.h"

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mi_lcao(const int& step, bool print)
{
    lcao::cal_mi_lcao(this->state_, this->p_operator, this->dm_, step, print);
}

// cal_mi_lcao<double> stub lives in template_helpers.cpp (single definition).

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::set_operator(
    hamilt::Operator<std::complex<double>>* op_in)
{
    this->p_operator = op_in;
}

template <>
void spinconstrain::SpinConstrain<double>::set_operator(
    hamilt::Operator<double>* op_in)
{
    this->p_operator = op_in;
}

#endif // __LCAO
