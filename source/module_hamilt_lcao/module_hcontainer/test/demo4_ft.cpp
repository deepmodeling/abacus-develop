#include "../hcontainer.h"

void folding_HR(const hamilt::HContainer<double>& hR,
                std::complex<double>* hk,
                const ModuleBase::Vector3<double>& kvec_d_in,
                const int ncol)
{
    for (int i = 0; i < hR.size_R_loop(); ++i)
    {
        // get R index
        int rx, ry, rz;
        hR.loop_R(i, rx, ry, rz);
        // only deal with current_R for hR
        hR.fix_R(rx, ry, rz);

        // cal k_phase
        ModuleBase::Vector3<double> dR(rx, ry, rz);
        const double arg = (kvec_d_in * dR) * ModuleBase::TWO_PI;
        double sinp, cosp;
        ModuleBase::libm::sincos(arg, &sinp, &cosp);
        const std::complex<double> kphase = std::complex<double>(cosp, sinp);

        // loop_atom_pairs
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int j = 0; j < hR.size_atom_pair(); ++j)
        {
            // Hk += HR * e^ikR
            hR.get_atom_pair(j).add_to_matrix(hk, ncol, kphase, 0);
        }
    }
}