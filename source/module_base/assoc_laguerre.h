#ifndef ASSOC_LAGUEERRE_H
#define ASSOC_LAGUEERRE_H

#include <vector>
#include <string>

class Assoc_Laguerre
{
    public:
        Assoc_Laguerre();
        ~Assoc_Laguerre();
        /// @brief generate the associated Laguerre polynomial (overloaded for double*)
        /// @param n principal quantum number
        /// @param l orbital quantum number
        /// @param ns number of x-coordinates
        /// @param s x-coordinates
        /// @param L y-coordinates
        void generate(const int &n, const int &l, const double ns, double* const &s, double* L);
        /// @brief generate the associated Laguerre polynomial (overloaded for std::vector)
        /// @param n principal quantum number
        /// @param l orbital quantum number
        /// @param x x-coordinates in std::vector
        /// @param y y-coordinates in std::vector
        void generate(const int &n, const int &l, std::vector<double> &x, std::vector<double> &y);
        /// @brief Laguerre polynomial
        /// @param n degree of the polynomial
        /// @param x radial coordinate
        /// @return L_n(x)
        double laguerre(const int &n, const double x);
        /// @brief recursive relationship to find the associated Laguerre polynomial
        /// @param n degree of the polynomial
        /// @param x radial coordinate
        /// @param a order of the polynomial
        /// @return L^(a)_n(x)
        double associate_laguerre(const int &n, const double x, const int &a);
        /// @brief wrapper for associate_laguerre
        /// @param n principal quantum number
        /// @param l orbital quantum number
        /// @param s radial coordinate 
        /// @return L^(2l+1)_(n-l-1)(s)
        double value(const int &n, const int &l, const double &s);
        /// @brief factorial function
        /// @param n 
        /// @return n!
        int factorial(const int &n);
};
#endif // ASSOC_LAGUEERRE_H