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

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::calculate_MW(
    const std::vector<std::vector<std::vector<double>>>& AorbMulP)
{
    size_t nw = this->get_nw();
    int nat = this->get_nat();

    this->zero_Mi();

    const int nlocal = nw / 2;
    for (const auto& sc_elem: this->get_atomCounts())
    {
        int it = sc_elem.first;
        int nat_it = sc_elem.second;
        for (int ia = 0; ia < nat_it; ia++)
        {
            int num = 0;
            int iat = this->get_iat(it, ia);
            std::vector<double> total_charge_soc(this->nspin_, 0.0);
            for (const auto& lnchi: this->get_lnchiCounts().at(it))
            {
                std::vector<double> sum_l(this->nspin_, 0.0);
                int L = lnchi.first;
                int nchi = lnchi.second;
                for (int Z = 0; Z < nchi; ++Z)
                {
                    std::vector<double> sum_m(this->nspin_, 0.0);
                    for (int M = 0; M < (2 * L + 1); ++M)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            sum_m[j] += AorbMulP[j][iat][num];
                        }
                        num++;
                    }
                    for (int j = 0; j < 4; j++)
                    {
                        sum_l[j] += sum_m[j];
                    }
                }
                for (int j = 0; j < 4; j++)
                {
                    total_charge_soc[j] += sum_l[j];
                }
            }
            this->Mi_[iat].x = total_charge_soc[1];
            this->Mi_[iat].y = total_charge_soc[2];
            this->Mi_[iat].z = total_charge_soc[3];
            if (std::abs(this->Mi_[iat].x) < this->sc_thr_)
                this->Mi_[iat].x = 0.0;
            if (std::abs(this->Mi_[iat].y) < this->sc_thr_)
                this->Mi_[iat].y = 0.0;
            if (std::abs(this->Mi_[iat].z) < this->sc_thr_)
                this->Mi_[iat].z = 0.0;
        }
    }
}