/**
 * @file k_vector_utils.cpp
 * @brief Implementation of k-vector utility functions.
 * @author rhx (created on 25-6-3)
 *
 * @note Since 2026-08-14 these free functions are thin wrappers around the
 *       spin-free members of ModuleCell::ReciprocalGrid / the K_Vectors
 *       IBZ orchestration, so that existing call sites (esolver_fp.cpp,
 *       klist.cpp, tests) keep working unchanged.
 */
#include "k_vector_utils.h"

#include "klist.h"
#include "source_base/global_variable.h"
#include "source_base/matrix3.h"

#include "source_base/formatter.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"

namespace KVectorUtils
{
void kvec_d2c(K_Vectors& kv, const ModuleBase::Matrix3& reciprocal_vec)
{
    kv.kvec_d2c(reciprocal_vec);
}
void kvec_c2d(K_Vectors& kv, const ModuleBase::Matrix3& latvec)
{
    kv.kvec_c2d(latvec);
}

void set_both_kvec(K_Vectors& kv, const ModuleBase::Matrix3& G, const ModuleBase::Matrix3& R, std::string& skpt)
{
    kv.set_both_kvec(G, R, skpt);
}

void set_after_vc(K_Vectors& kv, const int& nspin_in, const ModuleBase::Matrix3& reciprocal_vec)
{
    GlobalV::ofs_running << "\n SETUP K-POINTS" << std::endl;
    kv.set_nspin(nspin_in);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nspin", kv.get_nspin());

    // set cartesian k vectors.
    kv.kvec_d2c(reciprocal_vec);

    std::string table;
    table += "K-POINTS DIRECT COORDINATES\n";
    table += FmtCore::format("%8s%12s%12s%12s%8s\n", "KPOINTS", "DIRECT_X", "DIRECT_Y", "DIRECT_Z", "WEIGHT");
    for (int i = 0; i < kv.get_nks(); i++)
    {
        table += FmtCore::format("%8d%12.8f%12.8f%12.8f%8.4f\n",
                                 i + 1,
                                 kv.kvec_d[i].x,
                                 kv.kvec_d[i].y,
                                 kv.kvec_d[i].z,
                                 kv.wk[i]);
    }
    GlobalV::ofs_running << table << std::endl;

    kv.kd_done = true;
    kv.kc_done = true;

    print_klists(kv, GlobalV::ofs_running);
}

void print_klists(const K_Vectors& kv, std::ofstream& ofs)
{
    kv.print_klists(ofs);
}

#ifdef __MPI
void kvec_mpi_k(K_Vectors& kv)
{
    ModuleBase::TITLE("KVectorUtils", "kvec_mpi_k");

    Parallel_Common::bcast_bool(kv.kc_done);

    Parallel_Common::bcast_bool(kv.kd_done);

    Parallel_Common::bcast_int(kv.nspin);

    Parallel_Common::bcast_int(kv.nkstot);

    Parallel_Common::bcast_int(kv.nkstot_full);

    Parallel_Common::bcast_int(kv.nmp, 3);

    kv.kl_segids.resize(kv.nkstot);
    Parallel_Common::bcast_int(kv.kl_segids.data(), kv.nkstot);

    Parallel_Common::bcast_double(kv.koffset, 3);

    kv.nks = kv.para_k.nks_pool[GlobalV::MY_POOL];

    GlobalV::ofs_running << std::endl;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Number of k-points in this process", kv.nks);
    int nks_minimum = kv.nks;

    Parallel_Reduce::reduce_min(nks_minimum);

    if (nks_minimum == 0)
    {
        ModuleBase::WARNING_QUIT("K_Vectors::mpi_k()", " nks == 0, some processor have no k points!");
    }
    else
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Minimum distributed k-point number", nks_minimum);
    }

    std::vector<int> isk_aux(kv.nkstot);
    std::vector<double> wk_aux(kv.nkstot);
    std::vector<double> kvec_c_aux(kv.nkstot * 3);
    std::vector<double> kvec_d_aux(kv.nkstot * 3);
    std::vector<double> kvec_c_full_aux(kv.nkstot_full * 3);

    // collect and process in rank 0
    if (GlobalV::MY_RANK == 0)
    {
        for (int ik = 0; ik < kv.nkstot; ik++)
        {
            isk_aux[ik] = kv.isk[ik];
            wk_aux[ik] = kv.wk[ik];
            kvec_c_aux[3 * ik] = kv.kvec_c[ik].x;
            kvec_c_aux[3 * ik + 1] = kv.kvec_c[ik].y;
            kvec_c_aux[3 * ik + 2] = kv.kvec_c[ik].z;
            kvec_d_aux[3 * ik] = kv.kvec_d[ik].x;
            kvec_d_aux[3 * ik + 1] = kv.kvec_d[ik].y;
            kvec_d_aux[3 * ik + 2] = kv.kvec_d[ik].z;
            kvec_c_full_aux[3 * ik] = kv.kvec_c_full[ik].x;
            kvec_c_full_aux[3 * ik + 1] = kv.kvec_c_full[ik].y;
            kvec_c_full_aux[3 * ik + 2] = kv.kvec_c_full[ik].z;
        }
    }

    // broadcast k point data to all processors
    Parallel_Common::bcast_int(isk_aux.data(), kv.nkstot);

    Parallel_Common::bcast_double(wk_aux.data(), kv.nkstot);
    Parallel_Common::bcast_double(kvec_c_aux.data(), kv.nkstot * 3);
    Parallel_Common::bcast_double(kvec_d_aux.data(), kv.nkstot * 3);
    Parallel_Common::bcast_double(kvec_c_full_aux.data(), kv.nkstot_full * 3);

    // process k point data in each processor
    kv.renew(kv.nks * kv.nspin);

    // distribute
    int k_index = 0;

    for (int i = 0; i < kv.nks; i++)
    {
        // 3 is because each k point has three value:kx, ky, kz
        k_index = i + kv.para_k.startk_pool[GlobalV::MY_POOL];
        kv.kvec_c[i].x = kvec_c_aux[k_index * 3];
        kv.kvec_c[i].y = kvec_c_aux[k_index * 3 + 1];
        kv.kvec_c[i].z = kvec_c_aux[k_index * 3 + 2];
        kv.kvec_d[i].x = kvec_d_aux[k_index * 3];
        kv.kvec_d[i].y = kvec_d_aux[k_index * 3 + 1];
        kv.kvec_d[i].z = kvec_d_aux[k_index * 3 + 2];
        kv.kvec_c_full[i].x = kvec_c_full_aux[k_index * 3];
        kv.kvec_c_full[i].y = kvec_c_full_aux[k_index * 3 + 1];
        kv.kvec_c_full[i].z = kvec_c_full_aux[k_index * 3 + 2];
        kv.wk[i] = wk_aux[k_index];
        kv.isk[i] = isk_aux[k_index];
    }

#ifdef __EXX
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    { // bcast kstars
        kv.kstars.resize(kv.nkstot);
        for (int ikibz = 0; ikibz < kv.nkstot; ++ikibz)
        {
            int starsize = kv.kstars[ikibz].size();
            Parallel_Common::bcast_int(starsize);
            auto ks = kv.kstars[ikibz].begin();
            for (int ik = 0; ik < starsize; ++ik)
            {
                int isym = 0;
                ModuleBase::Vector3<double> ks_vec(0, 0, 0);
                if (GlobalV::MY_RANK == 0)
                {
                    isym = ks->first;
                    ks_vec = ks->second;
                    ++ks;
                }
                Parallel_Common::bcast_int(isym);
                Parallel_Common::bcast_double(ks_vec.x);
                Parallel_Common::bcast_double(ks_vec.y);
                Parallel_Common::bcast_double(ks_vec.z);
                if (GlobalV::MY_RANK != 0)
                {
                    kv.kstars[ikibz].insert(std::make_pair(isym, ks_vec));
                }
            }
        }
    }
#endif
} // END SUBROUTINE
#endif

void kvec_ibz_kpoint(K_Vectors& kv,
                     const ModuleSymmetry::Symmetry& symm,
                     bool use_symm,
                     std::string& skpt,
                     const UnitCell& ucell,
                     bool& match)
{
    kv.reduce_by_symmetry(ucell, symm, use_symm, skpt, match);
}
} // namespace KVectorUtils
