#include "to_w90_lcao.h"

#include "source_lcao/module_ri/abfs_vector3_order.h"
#include "source_io/module_parameter/parameter.h"
#include "fr_overlap.h"
#include "source_base/math_integral.h"
#include "source_base/math_polyint.h"
#include "source_base/math_sphbes.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_hamilt/module_hcontainer/atom_pair.h"

#include <fstream>
#include <functional>

#ifdef __LCAO

void toW90_LCAO::cal_Mmn(const UnitCell& ucell, const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi)
{
    // write .mmn file
    std::ofstream mmn_file;

    if (GlobalV::MY_RANK == 0)
    {
        std::string fileaddress = PARAM.globalv.global_out_dir + wannier_file_name + ".mmn";
        mmn_file.open(fileaddress.c_str(), std::ios::out);

        time_t time_now = time(nullptr);
        mmn_file << " Created on " << ctime(&time_now);
        mmn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << nntot << std::endl;
    }

    for (int ik = 0; ik < cal_num_kpts; ik++)
    {
        for (int ib = 0; ib < nntot; ib++)
        {
            int ikb = nnlist[ik][ib];
            ModuleBase::Vector3<double> phase_G = nncell[ik][ib];
            ModuleBase::ComplexMatrix Mmn;

            int cal_ik = ik + start_k_index;
            int cal_ikb = ikb + start_k_index;
            unkdotkb(ucell,kv, psi, cal_ik, cal_ikb, phase_G, Mmn);

            if (GlobalV::MY_RANK == 0)
            {
                mmn_file << std::setw(5) << ik + 1 << std::setw(5) << ikb + 1 << std::setw(5) << int(phase_G.x)
                         << std::setw(5) << int(phase_G.y) << std::setw(5) << int(phase_G.z) << std::endl;

                for (int n = 0; n < num_bands; n++)
                {
                    for (int m = 0; m < num_bands; m++)
                    {
                        mmn_file << std::setw(18) << std::setprecision(12) << std::showpoint << std::fixed
                                 << Mmn(m, n).real() << std::setw(18) << std::setprecision(12) << std::showpoint
                                 << std::fixed
                                 << Mmn(m, n).imag()
                                 // jingan test
                                 // << "    " << std::setw(12) << std::setprecision(9) << std::abs(Mmn(m, n))
                                 << std::endl;
                    }
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0) {
        mmn_file.close();
}
}

void toW90_LCAO::cal_Amn(const UnitCell& ucell, const K_Vectors& kv, const psi::Psi<std::complex<double>>& psi)
{
    produce_trial_in_lcao();
    construct_overlap_table_project();
    cal_orbA_overlap_R(ucell);

    // write .amn file
    std::ofstream Amn_file;

    if (GlobalV::MY_RANK == 0)
    {
        time_t time_now = time(nullptr);
        std::string fileaddress = PARAM.globalv.global_out_dir + wannier_file_name + ".amn";
        Amn_file.open(fileaddress.c_str(), std::ios::out);
        Amn_file << " Created on " << ctime(&time_now);
        Amn_file << std::setw(12) << num_bands << std::setw(12) << cal_num_kpts << std::setw(12) << num_wannier
                 << std::endl;
    }

    for (int ik = start_k_index; ik < (cal_num_kpts + start_k_index); ik++)
    {
        ModuleBase::ComplexMatrix Amn;
        unkdotA(kv, psi, ik, Amn);

        if (GlobalV::MY_RANK == 0)
        {
            for (int iw = 0; iw < num_wannier; iw++)
            {
                for (int ib = 0; ib < num_bands; ib++)
                {
                    Amn_file << std::setw(5) << ib + 1 << std::setw(5) << iw + 1 << std::setw(5)
                             << ik + 1 - start_k_index << std::setw(18) << std::showpoint << std::fixed
                             << std::setprecision(12) << Amn(ib, iw).real() << std::setw(18) << std::showpoint
                             << std::fixed << std::setprecision(12)
                             << Amn(ib, iw).imag()
                             // jingan test
                             //<< "   " << std::setw(18) << std::setprecision(13) << std::abs(Amn(ib, iw))
                             << std::endl;
                }
            }
        }
    }

    if (GlobalV::MY_RANK == 0) {
        Amn_file.close();
}
}

void toW90_LCAO::out_unk(const psi::Psi<std::complex<double>>& psi)
{
}

#endif
