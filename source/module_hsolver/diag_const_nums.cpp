#include "diag_const_nums.h"

#include <complex>

template class const_nums<double>;
template class const_nums<float>;
template class const_nums<std::complex<double>>;
template class const_nums<std::complex<float>>;

// Specialize templates to support double types
template <>
const_nums<double>::const_nums() : zero(0.0), one(1.0), neg_one(-1.0)
{
}

// Specialize templates to support double types
template <>
const_nums<float>::const_nums() : zero(0.0), one(1.0), neg_one(-1.0)
{
}

// Specialized templates to support std:: complex<double>types
template <>
const_nums<std::complex<double>>::const_nums()
    : zero(std::complex<double>(0.0, 0.0)), one(std::complex<double>(1.0, 0.0)),
      neg_one(std::complex<double>(-1.0, 0.0))
{
}

// Specialized templates to support std:: complex<float>types
template <>
const_nums<std::complex<float>>::const_nums()
    : zero(std::complex<float>(0.0, 0.0)), one(std::complex<float>(1.0, 0.0)), neg_one(std::complex<float>(-1.0, 0.0))
{
}
