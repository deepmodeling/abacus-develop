#include "spin_constrain.h"

template <>
std::vector<std::vector<std::vector<double>>> SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::convert(
    const ModuleBase::matrix& orbMulP)
{
    std::vector<std::vector<std::vector<double>>> AorbMulP;
    AorbMulP.resize(this->nspin_);
    int nat = this->get_nat();
    for (int is = 0; is < this->nspin_; ++is)
    {
        int num = 0;
        AorbMulP[is].resize(nat);
        for (const auto& sc_elem: this->get_atomCounts())
        {
            int it = sc_elem.first;
            int nat_it = sc_elem.second;
            int nw_it = this->get_orbitalCounts().at(it);
            for (int ia = 0; ia < nat_it; ia++)
            {
                int iat = this->get_iat(it, ia);
                AorbMulP[is][iat].resize(nw_it, 0.0);
                for (int iw = 0; iw < nw_it; iw++)
                {
                    AorbMulP[is][iat][iw] = orbMulP(is, num);
                    num++;
                }
            }
        }
    }

    return AorbMulP;
}