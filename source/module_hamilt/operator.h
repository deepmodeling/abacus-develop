#ifndef __OPERATOR
#define __OPERATOR

#include<complex>

namespace hamilt
{

class Operator
{
    public:
    virtual void act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size)const {return;}
    virtual void act(std::complex<double> *hk_matrix)const {return;}
    virtual void act(double *hk_matrix)const {return;}

    virtual void init(const int ik){this->ik = ik;}

    protected:
    int ik = 0;
};

}//end namespace hamilt

#endif