#ifdef __EXX

#include "spar_exx.h"

// --------------------------------------------------------
// Header files needed for implementation only
// --------------------------------------------------------

#include <RI/global/Global_Func-2.h>
#include <RI/ri/Cell_Nearest.h>
#include <stdexcept>

#include "source_base/abfs-vector3_order.h"
#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_ri/RI_2D_Comm.h"
#include "source_lcao/module_ri/RI_Util.hpp"

// --------------------------------------------------------
// Implementation of the cal_HR_exx function
// --------------------------------------------------------

namespace sparse_format
{

/**
 * @brief Implementation of the cal_HR_exx function
 */
template <typename Tdata>
void cal_HR_exx(
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const int& current_spin,
    const double& sparse_threshold,
    const int (&nmp)[3],
    const std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<Tdata>>>>& Hexxs)
{
    ModuleBase::TITLE("sparse_format", "cal_HR_exx");
    ModuleBase::timer::tick("sparse_format", "cal_HR_exx");

    const Tdata frac = GlobalC::exx_info.info_global.hybrid_alpha;

    std::map<int, std::array<double, 3>> atoms_pos;
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        atoms_pos[iat] = RI_Util::Vector3_to_array3(ucell.atoms[ucell.iat2it[iat]].tau[ucell.iat2ia[iat]]);
    }
    const std::array<std::array<double, 3>, 3> latvec
        = {RI_Util::Vector3_to_array3(ucell.a1),
           RI_Util::Vector3_to_array3(ucell.a2),
           RI_Util::Vector3_to_array3(ucell.a3)};

    const std::array<int, 3> Rs_period = {nmp[0], nmp[1], nmp[2]};

    RI::Cell_Nearest<int, int, 3, double, 3> cell_nearest;
    cell_nearest.init(atoms_pos, latvec, Rs_period);

    const std::vector<int> is_list
        = (PARAM.inp.nspin != 4) ? std::vector<int>{current_spin} : std::vector<int>{0, 1, 2, 3};

    for (const int is: is_list)
    {
        int is0_b = 0;
        int is1_b = 0;
        std::tie(is0_b, is1_b) = RI_2D_Comm::split_is_block(is);

        if (Hexxs.empty())
        {
            break;
        }

        for (const auto& HexxA: Hexxs[is])
        {
            const int iat0 = HexxA.first;
            for (const auto& HexxB: HexxA.second)
            {
                const int iat1 = HexxB.first.first;

                const Abfs::Vector3_Order<int> R = RI_Util::array3_to_Vector3(
                    cell_nearest.get_cell_nearest_discrete(iat0, iat1, HexxB.first.second));

                HS_Arrays.all_R_coor.insert(R);

                const RI::Tensor<Tdata>& Hexx = HexxB.second;

                for (size_t iw0 = 0; iw0 < Hexx.shape[0]; ++iw0)
                {
                    const int iwt0 = RI_2D_Comm::get_iwt(ucell, iat0, iw0, is0_b);
                    const int iwt0_local = pv.global2local_row(iwt0);

                    if (iwt0_local < 0)
                    {
                        continue;
                    }

                    for (size_t iw1 = 0; iw1 < Hexx.shape[1]; ++iw1)
                    {
                        const int iwt1 = RI_2D_Comm::get_iwt(ucell, iat1, iw1, is1_b);
                        const int iwt1_local = pv.global2local_col(iwt1);

                        if (iwt1_local < 0)
                        {
                            continue;
                        }

                        if (std::abs(Hexx(iw0, iw1)) > sparse_threshold)
                        {
                            if (PARAM.inp.nspin == 1 || PARAM.inp.nspin == 2)
                            {
                                auto& HR_sparse_ptr = HS_Arrays.HR_sparse[current_spin][R][iwt0];
                                double& HR_sparse = HR_sparse_ptr[iwt1];
                                HR_sparse += RI::Global_Func::convert<double>(frac * Hexx(iw0, iw1));
                                if (std::abs(HR_sparse) <= sparse_threshold)
                                {
                                    HR_sparse_ptr.erase(iwt1);
                                }
                            }
                            else if (PARAM.inp.nspin == 4)
                            {
                                auto& HR_sparse_ptr = HS_Arrays.HR_soc_sparse[R][iwt0];

                                std::complex<double>& HR_sparse = HR_sparse_ptr[iwt1];

                                HR_sparse += RI::Global_Func::convert<std::complex<double>>(frac * Hexx(iw0, iw1));

                                if (std::abs(HR_sparse) <= sparse_threshold)
                                {
                                    HR_sparse_ptr.erase(iwt1);
                                }
                            }
                            else
                            {
                                throw std::invalid_argument(std::string(__FILE__) + " line "
                                                            + std::to_string(__LINE__));
                            }
                        }
                    }
                }
            }
        }
    }

    ModuleBase::timer::tick("sparse_format", "cal_HR_exx");
}

// --------------------------------------------------------
// Explicit instantiations for double and complex<double> types
// --------------------------------------------------------

template void cal_HR_exx<double>(
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const int& current_spin,
    const double& sparse_thr,
    const int (&nmp)[3],
    const std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>>& Hexxs);

template void cal_HR_exx<std::complex<double>>(
    const UnitCell& ucell,
    const Parallel_Orbitals& pv,
    LCAO_HS_Arrays& HS_Arrays,
    const int& current_spin,
    const double& sparse_thr,
    const int (&nmp)[3],
    const std::vector<std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<std::complex<double>>>>>& Hexxs);

} // namespace sparse_format

#endif // __EXX

