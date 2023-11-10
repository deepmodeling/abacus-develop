#include "diagh.h"
template class consts<double>;
template class consts<std::complex<double>>;
template class consts<std::complex<float>>;

template<> consts<double>::consts()
{
    zero = 0.0;
    one = 1.0;
    neg_one = -1.0;
}
template<> consts<std::complex<double>>::consts()
{
    zero = std::complex<double>(0.0, 0.0);
    one = std::complex<double>(1.0, 0.0);
    neg_one = std::complex<double>(-1.0, 0.0);
}
template<> consts<std::complex<float>>::consts()
{
    zero = std::complex<float>(0.0, 0.0);
    one = std::complex<float>(1.0, 0.0);
    neg_one = std::complex<float>(-1.0, 0.0);
}
