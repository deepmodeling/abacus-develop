#include "module_base/assoc_laguerre.h"
#include "module_base/global_function.h"
//#include <tr1/cmath> // use cmath the factorial function
#include <cmath>
Assoc_Laguerre::Assoc_Laguerre()
{
}

Assoc_Laguerre::~Assoc_Laguerre()
{
}

void Assoc_Laguerre::generate(const int &n, const int &l, const double ns, double* const &s, double* L)
{
    for(int i = 0; i < ns; i++)
    {
        L[i] = this->value(n, l, s[i]);
    }
}

void Assoc_Laguerre::generate(const int &n, const int &l, std::vector<double> &x, std::vector<double> &y)
{
    for(int i = 0; i < x.size(); i++)
    {
        y[i] = this->value(n, l, x[i]);
    }
}

double Assoc_Laguerre::laguerre(const int &n, const double x)
{
    if(n == 0)
    {
        return 1;
    }
    else if(n == 1)
    {
        return -x + 1;
    }
    else if(n == 2)
    {
        return 0.5 * x * x - 2 * x + 1;
    }
    else if(n == 3)
    {
        return -x * x * x / 6.0 + 3.0 * x * x / 2.0 - 3.0 * x + 1;
    }
    else if(n >= 4)
    {
        double n_ = static_cast<double>(n);
        double first = (2*n_ - 1 - x)/n_ * Assoc_Laguerre::laguerre(n-1, x);
        double second = (n_ - 1)/n_ * Assoc_Laguerre::laguerre(n-2, x);
        return first - second;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Assoc_Laguerre::laguerre", "n is out of range");
        return 0;
    }
}

double Assoc_Laguerre::associate_laguerre(const int &n, const double x, const int &a)
{
    // formula from https://en.wikipedia.org/wiki/Laguerre_polynomials
    double n_ = static_cast<double>(n);
    double a_ = static_cast<double>(a);
    if(n == 0)
    {
        return 1;
    }
    else if(n == 1)
    {
        return -x + 1 + a_;
    }
    else if(n == 2)
    {
        return 0.5 * (x*x - 2*(a_+2)*x + (a_+1)*(a_+2));
    }
    else if(n == 3)
    {
        return -x*x*x/6.0 + (a_+3)*x*x/2.0 - (a_+2)*(a_+3)*x/2.0 + (a_+1)*(a_+2)*(a_+3)/6.0;
    }
    else if(n >= 4)
    {
        double first = (2*n_ - 1 + a_ - x)/n_ * this->associate_laguerre(n-1, x, a);
        double second = (n_ + a_ - 1)/n_ * this->associate_laguerre(n-2, x, a);
        return first - second;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Assoc_Laguerre::associate_laguerre", "n is out of range");
        return 0;
    }
}

int Assoc_Laguerre::factorial(const int &n)
{
    if(n == 0)
    {
        return 1;
    }
    else if(n > 0)
    {
        return n * this->factorial(n-1);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Assoc_Laguerre::factorial", "n is out of range");
        return 0;
    }
}

double Assoc_Laguerre::value(const int &n, const int &l, const double &s)
{
    int k_ = 2*l + 1;
    int n_ = n - l - 1;
    if(k_ < 0)
    {
        ModuleBase::WARNING_QUIT("Assoc_Laguerre::value", "k is out of range");
        return 0;
    }
    if(n_ < 0)
    {
        ModuleBase::WARNING_QUIT("Assoc_Laguerre::value", "n is out of range");
        return 0;
    }
    double L = 0;
    for(int iq = 0; iq <= n_; iq++)
    {
        L += std::pow(-s, iq) * 
        static_cast<double>(this->factorial(n_ + k_)) / 
        static_cast<double>(this->factorial(n_ - iq)) / 
        static_cast<double>(this->factorial(k_ + iq)) / 
        static_cast<double>(this->factorial(iq));
    }
    //L = std::tr1::assoc_laguerre(n_, k_, s); // use standard library
    return L;
}