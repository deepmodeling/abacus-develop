#ifndef FUNC_HPP
#define FUNC_HPP
template <class T, class TI>
inline void ZEROS(std::complex<T>* u, const TI n) // Peize Lin change int to TI at 2020.03.03
{
    assert(n >= 0);
    for (TI i = 0; i < n; i++)
    {
        u[i] = std::complex<T>(0.0, 0.0);
    }
    return;
}

template <class T, class TI>
inline void ZEROS(T* u, const TI n) // Peize Lin change int to TI at 2020.03.03
{
    assert(n >= 0);
    for (TI i = 0; i < n; i++)
    {
        u[i] = 0;
    }
}
#endif