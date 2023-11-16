#include "associate_lagurre.h"

template <typename T>
Associate_Laguerre<T>::Associate_Laguerre()
{
}

template <typename T>
Associate_Laguerre<T>::~Associate_Laguerre()
{
}

template <typename T>
T Associate_Laguerre<T>::value(const int& n, const int& alpha, T x)
{
    if(n < 0)
    {
        return 0;
    }
    else if(n == 0)
    {
        return 1;
    }
    else if(n == 1)
    {
        return 1 + alpha - x;
    }
    else
    {
        return ((2*n + alpha - 1 - x) * Associate_Laguerre::value(n - 1, alpha, x) - (n + alpha - 1) * Associate_Laguerre::value(n - 2, alpha, x)) / n;
    }
}

template <typename T>
void Associate_Laguerre<T>::generate(const int& n, const int& alpha, T* x, T* Laguerre, const int& size)
{
    for(int i = 0; i < size; i++)
    {
        Laguerre[i] = Associate_Laguerre::value(n, alpha, x[i]);
    }
}

