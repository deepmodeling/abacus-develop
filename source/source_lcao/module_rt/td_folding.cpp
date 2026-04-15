#include "td_folding.h"
#include "source_base/libm/libm.h"
namespace module_rt{
template<typename TR>
void folding_HR_td(const UnitCell& ucell,
                const hamilt::HContainer<TR>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const ModuleBase::Vector3<double>& cart_At,
                const int ncol,
                const int hk_type)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < hR.size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR>& tmp = hR.get_atom_pair(i);
        for(int ir = 0;ir < tmp.get_R_size(); ++ir )
        {
            const ModuleBase::Vector3<int> r_index = tmp.get_R_index(ir);

            //new
            //cal tddft phase for hybrid gauge
            const int iat1 = tmp.get_atom_i();
            const int iat2 = tmp.get_atom_j();
            ModuleBase::Vector3<double> dtau = ucell.cal_dtau(iat1, iat2, r_index);
            const double arg_td = cart_At * dtau * ucell.lat0;
            //new

            // cal k_phase
            // if TK==std::complex<double>, kphase is e^{ikR}
            const ModuleBase::Vector3<double> dR(r_index.x, r_index.y, r_index.z);
            const double arg = (kvec_d_in * dR) * ModuleBase::TWO_PI + arg_td;
            double sinp = 0.0, cosp = 0.0;
            ModuleBase::libm::sincos(arg, &sinp, &cosp);
            std::complex<double> kphase = std::complex<double>(cosp, sinp);

            tmp.find_R(r_index);
            tmp.add_to_matrix(hk, ncol, kphase, hk_type);
        }
    }
}

template<typename TR>
void folding_partial_HR(const UnitCell& ucell,
                const hamilt::HContainer<TR>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ix,
                const int ncol,
                const int hk_type)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < hR.size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR>& tmp = hR.get_atom_pair(i);
        for(int ir = 0;ir < tmp.get_R_size(); ++ir )
        {
            const ModuleBase::Vector3<int> r_index = tmp.get_R_index(ir);
            const ModuleBase::Vector3<double> dR(r_index.x, r_index.y, r_index.z);
            const double arg = (kvec_d_in * dR) * ModuleBase::TWO_PI;
            double sinp = 0.0, cosp = 0.0;
            ModuleBase::libm::sincos(arg, &sinp, &cosp);
            std::complex<double> kphase = std::complex<double>(cosp, sinp);
            const ModuleBase::Vector3<double> dR_car = dR * ucell.latvec * ucell.lat0;

            tmp.find_R(r_index);
            tmp.add_to_matrix(hk, ncol, kphase * ModuleBase::IMAG_UNIT * std::complex<double>(dR_car[ix]), hk_type);
        }
    }
}

template<typename TR>
void folding_partial_HR_td(const UnitCell& ucell,
                        const hamilt::HContainer<TR>& hR,
                        std::complex<double>* hk,
                        const ModuleBase::Vector3<double>& kvec_d_in,
                        const ModuleBase::Vector3<double>& cart_At,
                        const int ix,
                        const int ncol,
                        const int hk_type)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < hR.size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<TR>& tmp = hR.get_atom_pair(i);
        for(int ir = 0;ir < tmp.get_R_size(); ++ir )
        {
            const ModuleBase::Vector3<int> r_index = tmp.get_R_index(ir);

            //new
            //cal tddft phase for mixing gague
            const int iat1 = tmp.get_atom_i();
            const int iat2 = tmp.get_atom_j();
            ModuleBase::Vector3<double> dtau = ucell.cal_dtau(iat1, iat2, r_index);
            const double arg_td = cart_At * dtau * ucell.lat0;

            //new
            // cal k_phase
            // if TK==std::complex<double>, kphase is e^{ikR}
            const ModuleBase::Vector3<double> dR(r_index.x, r_index.y, r_index.z);
            const double arg = (kvec_d_in * dR) * ModuleBase::TWO_PI + arg_td;
            double sinp = 0.0, cosp = 0.0;
            ModuleBase::libm::sincos(arg, &sinp, &cosp);
            std::complex<double> kphase = std::complex<double>(cosp, sinp);
            const ModuleBase::Vector3<double> dR_car = dR * ucell.latvec * ucell.lat0;

            tmp.find_R(r_index);
            tmp.add_to_matrix(hk, ncol, kphase * ModuleBase::IMAG_UNIT * std::complex<double>(dR_car[ix]), hk_type);
        }
    }
}
template
void folding_HR_td<double>(const UnitCell& ucell,
                const hamilt::HContainer<double>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const ModuleBase::Vector3<double>& At,
                const int ncol,
                const int hk_type);
template
void folding_HR_td<std::complex<double>>(const UnitCell& ucell,
                const hamilt::HContainer<std::complex<double>>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const ModuleBase::Vector3<double>& At,
                const int ncol,
                const int hk_type);
template
void folding_partial_HR<std::complex<double>>(const UnitCell& ucell,
                const hamilt::HContainer<std::complex<double>>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ix,
                const int ncol,
                const int hk_type);
template
void folding_partial_HR<double>(const UnitCell& ucell,
                const hamilt::HContainer<double>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ix,
                const int ncol,
                const int hk_type);
template
void folding_partial_HR_td<std::complex<double>>(const UnitCell& ucell,
                const hamilt::HContainer<std::complex<double>>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const ModuleBase::Vector3<double>& cart_At,
                const int ix,
                const int ncol,
                const int hk_type);
template
void folding_partial_HR_td<double>(const UnitCell& ucell,
                const hamilt::HContainer<double>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const ModuleBase::Vector3<double>& cart_At,
                const int ix,
                const int ncol,
                const int hk_type);
}// namespace module_rt