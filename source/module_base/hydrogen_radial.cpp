#include "hydrogen_radial.h"
#include <cmath>

template<typename T>
T Hydrogen_Radial<T>::value(const int& Z, const int& n, const int& l, const T& r)
{
    T a = hbar*hbar/(me*e*e);
    T rho = 2*Z*r/(n*a);
    T Rnl = std::sqrt(
        4*std::pow(Z, 3)/std::pow(n*a, 3)/n
        *(T)(factorial(n - l - 1)/factorial(n + l))
        )
        *
        std::exp(-rho/2.)
        *
        std::pow(rho, l)
        *
        L.value(n - l - 1, 2*l + 1, rho);

    return Rnl;
}

template<typename T>
void Hydrogen_Radial<T>::generate(const int& Z, const int& n, const int& l, T* r, T* Rnl, const int& size)
{
    for(int i = 0; i < size; i++)
    {
        Rnl[i] = Hydrogen_Radial::value(Z, n, l, r[i]);
    }
}