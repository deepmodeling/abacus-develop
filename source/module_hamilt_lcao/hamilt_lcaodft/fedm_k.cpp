#include "FORCE.h"
#include "module_base/memory.h"
#include "module_parameter/parameter.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_base/tool_threading.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_elecstate/cal_dm.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/module_dm/cal_dm_psi.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/write_HS.h"

#include <map>
#include <unordered_map>

#ifdef __DEEPKS
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"

template <>
void Force_LCAO<std::complex<double>>::cal_fedm(const bool isforce,
                                                const bool isstress,
                                                ForceStressArrays& fsr,
                                                const UnitCell& ucell,
                                                const elecstate::DensityMatrix<std::complex<double>, double>* dm,
                                                const psi::Psi<std::complex<double>>* psi,
                                                const Parallel_Orbitals& pv,
                                                const elecstate::ElecState* pelec,
                                                ModuleBase::matrix& foverlap,
                                                ModuleBase::matrix& soverlap,
                                                const K_Vectors* kv,
                                                Record_adj* ra)
{
    ModuleBase::TITLE("Force_LCAO", "cal_fedm");
    ModuleBase::timer::tick("Force_LCAO", "cal_fedm");

    const int nspin = PARAM.inp.nspin;
    const int nbands = PARAM.inp.nbands;

    // construct a DensityMatrix object
    elecstate::DensityMatrix<std::complex<double>, double> edm(kv, &pv, nspin);

    //--------------------------------------------
    // calculate the energy density matrix here.
    //--------------------------------------------

    ModuleBase::matrix wg_ekb;
    wg_ekb.create(kv->get_nks(), nbands);
    ModuleBase::Memory::record("Force::wg_ekb", sizeof(double) * kv->get_nks() * nbands);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 1024)
#endif
    for (int ik = 0; ik < kv->get_nks(); ik++)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            wg_ekb(ik, ib) = pelec->wg(ik, ib) * pelec->ekb(ik, ib);
        }
    }
    std::vector<ModuleBase::ComplexMatrix> edm_k(kv->get_nks());

    // use the original formula (Hamiltonian matrix) to calculate energy density matrix
    if (dm->EDMK.size())
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int ik = 0; ik < kv->get_nks(); ++ik)
        {
            edm.set_DMK_pointer(ik, dm->EDMK[ik].c);
        }
    }
    else
    {
        // cal_dm_psi
        elecstate::cal_dm_psi(edm.get_paraV_pointer(), wg_ekb, psi[0], edm);
    }

    // cal_dm_2d
    edm.init_DMR(*ra, &ucell);
    edm.cal_DMR();
    edm.sum_DMR_spin();
    //--------------------------------------------
    // summation \sum_{i,j} E(i,j)*dS(i,j)
    // BEGIN CALCULATION OF FORCE OF EACH ATOM
    //--------------------------------------------
    int total_irr = 0;
#ifdef _OPENMP
#pragma omp parallel
    {
        int num_threads = omp_get_num_threads();
        ModuleBase::matrix local_soverlap(3, 3);
        int local_total_irr = 0;
#else
    ModuleBase::matrix& local_soverlap = soverlap;
    int& local_total_irr = total_irr;
#endif

        ModuleBase::Vector3<double> tau1, dtau, tau2;

#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            const int T1 = ucell.iat2it[iat];
            Atom* atom1 = &ucell.atoms[T1];
            const int I1 = ucell.iat2ia[iat];
            // get iat1
            int iat1 = ucell.itia2iat(T1, I1);
            double* foverlap_iat;
            if (isforce)
            {
                foverlap_iat = &foverlap(iat, 0);
            }

#ifdef _OPENMP
            // using local stack to avoid false sharing in multi-threaded case
            double foverlap_temp[3] = {0.0, 0.0, 0.0};
            if (num_threads > 1)
            {
                foverlap_iat = foverlap_temp;
            }
#endif
            int irr = pv.nlocstart[iat];
            const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
            for (int cb = 0; cb < ra->na_each[iat]; ++cb)
            {
                const int T2 = ra->info[iat][cb][3];
                const int I2 = ra->info[iat][cb][4];

                const int start2 = ucell.itiaiw2iwt(T2, I2, 0);

                Atom* atom2 = &ucell.atoms[T2];

                // get iat2
                int iat2 = ucell.itia2iat(T2, I2);
                double Rx = ra->info[iat][cb][0];
                double Ry = ra->info[iat][cb][1];
                double Rz = ra->info[iat][cb][2];
                // get BaseMatrix
                hamilt::BaseMatrix<double>* tmp_matrix = edm.get_DMR_pointer(1)->find_matrix(iat1, iat2, Rx, Ry, Rz);
                if (tmp_matrix == nullptr)
                {
                    continue;
                }
                int row_ap = pv.atom_begin_row[iat1];
                int col_ap = pv.atom_begin_col[iat2];
                // get DMR
                for (int mu = 0; mu < pv.get_row_size(iat1); ++mu)
                {
                    for (int nu = 0; nu < pv.get_col_size(iat2); ++nu)
                    {
                        // here do not sum over spin due to edm.sum_DMR_spin();
                        double edm2d1 = tmp_matrix->get_value(mu, nu);
                        double edm2d2 = 2.0 * edm2d1;

                        if (isforce)
                        {
                            foverlap_iat[0] -= edm2d2 * fsr.DSloc_Rx[irr];
                            foverlap_iat[1] -= edm2d2 * fsr.DSloc_Ry[irr];
                            foverlap_iat[2] -= edm2d2 * fsr.DSloc_Rz[irr];
                        }
                        if (isstress)
                        {
                            for (int ipol = 0; ipol < 3; ipol++)
                            {
                                local_soverlap(0, ipol) += edm2d1 * fsr.DSloc_Rx[irr] * fsr.DH_r[irr * 3 + ipol];
                                if (ipol < 1)
                                {
                                    continue;
                                }
                                local_soverlap(1, ipol) += edm2d1 * fsr.DSloc_Ry[irr] * fsr.DH_r[irr * 3 + ipol];
                                if (ipol < 2)
                                {
                                    continue;
                                }
                                local_soverlap(2, ipol) += edm2d1 * fsr.DSloc_Rz[irr] * fsr.DH_r[irr * 3 + ipol];
                            }
                        }
                        //}
                        ++local_total_irr;
                        ++irr;
                    } // end kk
                }     // end jj
            }         // end cb
#ifdef _OPENMP
            if (isforce && num_threads > 1)
            {
                foverlap(iat, 0) += foverlap_iat[0];
                foverlap(iat, 1) += foverlap_iat[1];
                foverlap(iat, 2) += foverlap_iat[2];
            }
#endif
        } // end iat
#ifdef _OPENMP
#pragma omp critical(cal_foverlap_k_reduce)
        {
            total_irr += local_total_irr;
            if (isstress)
            {
                for (int ipol = 0; ipol < 3; ipol++)
                {
                    soverlap(0, ipol) += local_soverlap(0, ipol);
                    if (ipol < 1)
                    {
                        continue;
                    }
                    soverlap(1, ipol) += local_soverlap(1, ipol);
                    if (ipol < 2)
                    {
                        continue;
                    }
                    soverlap(2, ipol) += local_soverlap(2, ipol);
                }
            }
        }
    }
#endif

    if (isstress)
    {
        StressTools::stress_fill(ucell.lat0, ucell.omega, soverlap);
    }

    if (total_irr != pv.nnr)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong irr", total_irr);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong pv.nnr", pv.nnr);
        ModuleBase::WARNING_QUIT("Force_LCAO::fedm_k", "irr!=pv.nnr");
    }

    ModuleBase::timer::tick("Force_LCAO", "cal_fedm");
    return;
}
