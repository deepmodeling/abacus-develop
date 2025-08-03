#include "lr_util_hcontainer.h"
namespace LR_Util
{
    void get_DMR_real_imag_part(const elecstate::DensityMatrix<std::complex<double>, std::complex<double>>& DMR,
        elecstate::DensityMatrix<std::complex<double>, double>& DMR_real,
        const int& nat,
        const char& type)
    {
        assert(DMR.get_DMR_vector().size() == DMR_real.get_DMR_vector().size());
        bool get_imag = (type == 'I' || type == 'i');
        for (int is = 0;is < DMR.get_DMR_vector().size();++is)
        {
            auto dr = DMR.get_DMR_vector()[is]; //get_DMR_pointer() has bug when is=0
            auto dr_real = DMR_real.get_DMR_vector()[is];
            assert(dr != nullptr);
            assert(dr_real != nullptr);
            for (int ia = 0;ia < nat;ia++) {
                for (int ja = 0;ja < nat;ja++)
                {
                    auto ap = dr->find_pair(ia, ja);
                    auto ap_real = dr_real->find_pair(ia, ja);
                    for (int iR = 0;iR < ap->get_R_size();++iR)
                    {
                        // R index may be different between the two HContainers, find by R value instead of R-index
                        auto dR = ap->get_R_index(iR);
                        auto ptr = ap->get_HR_values(iR).get_pointer();
                        auto ptr_real = ap_real->get_HR_values(dR.x, dR.y, dR.z).get_pointer();
                        for (int i = 0;i < ap->get_size();++i) { ptr_real[i] = (get_imag ? ptr[i].imag() : ptr[i].real()); }
                    }
                }
            }
        }
    }

    void set_HR_real_imag_part(const hamilt::HContainer<double>& HR_real,
        hamilt::HContainer<std::complex<double>>& HR,
        const int& nat,
        const char& type)
    {
        bool get_imag = (type == 'I' || type == 'i');
        for (int ia = 0;ia < nat;ia++) {
            for (int ja = 0;ja < nat;ja++)
            {
                auto ap = HR.find_pair(ia, ja);
                auto ap_real = HR_real.find_pair(ia, ja);
                for (int iR = 0;iR < ap->get_R_size();++iR)
                {
                    // R index may be different between the two HContainers, find by R value instead of R-index
                    auto dR = ap->get_R_index(iR);
                    auto ptr = ap->get_HR_values(iR).get_pointer();
                    auto ptr_real = ap_real->get_HR_values(dR.x, dR.y, dR.z).get_pointer();
                    for (int i = 0;i < ap->get_size();++i) { get_imag ? ptr[i].imag(ptr_real[i]) : ptr[i].real(ptr_real[i]); }
                }
            }
        }
    }
}