#ifndef __EKINETICPW
#define __EKINETICPW

#include "operator_pw.h"

namespace hamilt
{

#ifndef __EKINETICTEMPLATE
#define __EKINETICTEMPLATE

template<class T> class Ekinetic : public T 
{};

#endif

template<>
class Ekinetic<OperatorPW> : public OperatorPW
{
    public:
    Ekinetic(
        double tpiba2_in, 
        const double* gk2_in,
        const int gk2_dim_in
    );

    virtual ~Ekinetic(){};

    virtual void act
    (
        const psi::Psi<std::complex<double>> *psi_in, 
        const int n_npwx, 
        const std::complex<double>* tmpsi_in, 
        std::complex<double>* tmhpsi
    )const override;

    private:

    mutable int max_npw = 0;

    mutable int npol = 0;

    double tpiba2 = 0.0;

    const double* gk2 = nullptr;

    int gk2_dim = 0;

};

} // namespace hamilt

#endif