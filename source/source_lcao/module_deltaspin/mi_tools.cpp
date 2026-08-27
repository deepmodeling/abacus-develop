#include "mi_tools.h"

#include "source_base/tool_quit.h"

namespace spinconstrain
{

void accumulate_Mi_from_becp(const std::complex<double>* becp,
                            int nkb,
                            int nbands,
                            int npol,
                            int spin_sign,
                            const double* wg_ik,
                            const int* nh_iat,
                            std::vector<ModuleBase::Vector3<double>>& mi)
{
    if (becp == nullptr)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "becp is nullptr");
    }
    if (wg_ik == nullptr)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "wg_ik is nullptr");
    }
    if (nh_iat == nullptr)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "nh_iat is nullptr");
    }
    if (nkb <= 0)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "nkb must be positive");
    }
    if (nbands <= 0)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "nbands must be positive");
    }
    if (npol != 1 && npol != 2)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "npol must be 1 or 2");
    }
    if (spin_sign != -1 && spin_sign != 1)
    {
        ModuleBase::WARNING_QUIT("accumulate_Mi_from_becp", "spin_sign must be -1 or 1");
    }

    const std::complex<double> zero(0.0, 0.0);
    if (npol == 2)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            const double weight = wg_ik[ib];
            int begin_ih = 0;
            for (int iat = 0; iat < static_cast<int>(mi.size()); iat++)
            {
                std::complex<double> occ[4] = {zero, zero, zero, zero};
                const int nh = nh_iat[iat];
                for (int ih = 0; ih < nh; ih++)
                {
                    const int index = ib * 2 * nkb + begin_ih + ih;
                    occ[0] += conj(becp[index]) * becp[index];
                    occ[1] += conj(becp[index]) * becp[index + nkb];
                    occ[2] += conj(becp[index + nkb]) * becp[index];
                    occ[3] += conj(becp[index + nkb]) * becp[index + nkb];
                }
                mi[iat] += pauli_to_moment(occ, weight);
                begin_ih += nh;
            }
        }
    }
    else // npol == 1
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            const double weight = wg_ik[ib];
            int begin_ih = 0;
            for (int iat = 0; iat < static_cast<int>(mi.size()); iat++)
            {
                double occ = 0.0;
                const int nh = nh_iat[iat];
                for (int ih = 0; ih < nh; ih++)
                {
                    const int index = ib * nkb + begin_ih + ih;
                    occ += (conj(becp[index]) * becp[index]).real();
                }
                mi[iat].z += weight * occ * spin_sign;
                begin_ih += nh;
            }
        }
    }
}

} // namespace spinconstrain
