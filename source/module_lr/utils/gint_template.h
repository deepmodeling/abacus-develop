#pragma once
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
namespace LR
{
    template <typename T> struct TGint;
    template <>
    struct TGint<double> {
        using type = Gint_Gamma;
    };
    template <>
    struct TGint<std::complex<double>> {
        using type = Gint_k;
    };
}