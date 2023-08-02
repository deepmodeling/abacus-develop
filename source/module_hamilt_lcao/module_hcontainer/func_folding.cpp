#include "hcontainer_funcs.h"
#include "module_base/libm/libm.h"

namespace hamilt
{
/**
 * @brief calculate the Hk matrix with specific k vector
 * @param hR the HContainer of <I,J,R> atom pairs
 * @param hk the data pointer of Hk matrix, the size of hk would be nrow * ncol
 * @param kvec_d_in the k vector in Direct coordinate
 * @param hk_ld the leading dimension number of hk, ncol for row-major, nrow for column-major
 * @param hk_type the data-type of hk, 0 is row-major, 1 is column-major
*/
template<typename TR>
void folding_HR(const hamilt::HContainer<TR>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ncol,
                const int hk_type)
{
    for (int i = 0; i < hR.size_R_loop(); ++i)
    {
        // get R index
        int rx, ry, rz;
        hR.loop_R(i, rx, ry, rz);
        // only deal with current_R for hR
        hR.fix_R(rx, ry, rz);

        // cal k_phase
        // if TK==std::complex<double>, kphase is e^{ikR}
        const ModuleBase::Vector3<double> dR(rx, ry, rz);
        const double arg = (kvec_d_in * dR) * ModuleBase::TWO_PI;
        double sinp, cosp;
        ModuleBase::libm::sincos(arg, &sinp, &cosp);
        std::complex<double> kphase = std::complex<double>(cosp, sinp);

        // loop_atom_pairs
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j < hR.size_atom_pairs(); ++j)
        {
            // Hk += HR * e^ikR
            hR.get_atom_pair(j).add_to_matrix(hk, ncol, kphase, hk_type);
        }
    }
}

// template instantiation
template void folding_HR<std::complex<double>>(const hamilt::HContainer<std::complex<double>>& hR,
                                            std::complex<double>* hk,
                                            const ModuleBase::Vector3<double>& kvec_d_in,
                                            const int ncol,
                                            const int hk_type);
template void folding_HR<double>(const hamilt::HContainer<double>& hR,
                                std::complex<double>* hk,
                                const ModuleBase::Vector3<double>& kvec_d_in,
                                const int ncol,
                                const int hk_type);
// special case for double
void folding_HR(const hamilt::HContainer<double>& hR,
                double* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ncol,
                const int hk_type)
{
    for (int i = 0; i < hR.size_R_loop(); ++i)
    {
        // get R index
        int rx, ry, rz;
        hR.loop_R(i, rx, ry, rz);
        // only deal with current_R for hR
        hR.fix_R(rx, ry, rz);

        // cal k_phase
        // if TK==double, kphase is 1.0
        double kphase = 1.0;

        // loop_atom_pairs
#ifdef _OPENMP
#pragma omp for schedule(static, 16)
#endif
        for (int j = 0; j < hR.size_atom_pairs(); ++j)
        {
            // Hk += HR * e^ikR
            hR.get_atom_pair(j).add_to_matrix(hk, ncol, kphase, hk_type);
        }
    }
}

} // namespace hamilt