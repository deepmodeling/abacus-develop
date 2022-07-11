#ifndef __EKINETICPW
#define __EKINETICPW

#include "operator.h"

namespace hamilt
{

class EkineticPW : public Operator
{
    public:
    EkineticPW(
        int max_npw_in,
        int npol_in,
        double tpiba2_in,
        const int* ngk_in, 
        const double* gk2_in
    );

    void act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const override;

    private:

    int max_npw = 0;

    int npol = 0;

    double tpiba2 = 0.0;

    const int* ngk = nullptr;

    const double* gk2 = nullptr;
};

} // namespace hamilt

#endif