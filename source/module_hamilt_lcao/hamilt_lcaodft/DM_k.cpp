#include "local_orbital_charge.h"
#include "module_base/blas_connector.h"
#include "module_base/libm/libm.h"
#include "module_base/memory.h"
#include "module_base/parallel_common.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#ifdef __MKL
#include <mkl_service.h>
#endif

void Local_Orbital_Charge::allocate_DM_k(const int& nks, const int& nnrg)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "allocate_k");

    this->nnrg_now = nnrg;
    // xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
    if (GlobalV::OUT_LEVEL != "m")
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nnrg_last", nnrg_last);
    if (GlobalV::OUT_LEVEL != "m")
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "nnrg_now", nnrg_now);

    if (this->init_DM_R)
    {
        assert(nnrg_last > 0);
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            delete[] DM_R[is];
        }
        delete[] DM_R;
        init_DM_R = false;
    }

    if (nnrg_now > 0)
    {
        this->DM_R = new double*[GlobalV::NSPIN];
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            this->DM_R[is] = new double[nnrg_now];
            ModuleBase::GlobalFunc::ZEROS(DM_R[is], nnrg_now);
        }
        this->nnrg_last = nnrg_now;
        this->init_DM_R = true;
        ModuleBase::Memory::record("LOC::DM_R", sizeof(double) * GlobalV::NSPIN * nnrg_now);
    }
    else if (nnrg_now == 0)
    {
        this->init_DM_R = false;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Local_Orbital_Charge::allocate_k", "check init_DM_R.");
    }

    // Peize Lin test 2019-01-16
    this->init_dm_2d(nks);

    return;
}

#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
inline void cal_DM_ATOM(const Grid_Technique& gt,
                        const std::complex<double> fac,
                        const Record_adj& RA,
                        const int ia1,
                        const int iw1_lo,
                        const int nw1,
                        const int gstart,
                        std::complex<double>*** wfc_k_grid,
                        std::complex<double>* WFC_PHASE,
                        std::complex<double>** DM_ATOM,
                        const ModuleBase::matrix& wg_in,
                        const K_Vectors& kv)
{

    const char transa = 'N';
    const char transb = 'T';
    const std::complex<double> alpha = 1;
    const std::complex<double> beta = 1;

    for (int ik = 0; ik < kv.get_nks(); ik++)
    {
        std::complex<double>** wfc = wfc_k_grid[ik];
        const int ispin = kv.isk[ik];
        int atom2start = 0;

        for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
        {
            std::complex<double>* DM = &DM_ATOM[ispin][atom2start];
            const int T2 = RA.info[ia1][ia2][3];
            const int I2 = RA.info[ia1][ia2][4];
            Atom* atom2 = &GlobalC::ucell.atoms[T2];
            const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
            const int iw2_lo = gt.trace_lo[start2];
            const int nw2 = atom2->nw;
            std::complex<double> exp_R = ModuleBase::libm::exp(fac
                                                               * (kv.kvec_d[ik].x * RA.info[ia1][ia2][0]
                                                                  + kv.kvec_d[ik].y * RA.info[ia1][ia2][1]
                                                                  + kv.kvec_d[ik].z * RA.info[ia1][ia2][2]));

            // ModuleBase::GlobalFunc::ZEROS(WFC_PHASE, GlobalV::NBANDS*nw1);
            int ibStart = 0;
            int nRow = 0;
            for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
            {
                const double wg_local = wg_in(ik, ib);
                if (wg_local > 0 || GlobalV::ocp == 1)
                {
                    if (nRow == 0)
                        ibStart = ib;
                    const int iline = nRow * nw1;
                    std::complex<double> phase = exp_R * wg_local;
                    for (int iw1 = 0; iw1 < nw1; ++iw1)
                    {
                        WFC_PHASE[iline + iw1] = phase * conj(wfc[ib][iw1_lo + iw1]);
                    }
                    ++nRow;
                }
                else
                {
                    break;
                }
            } // ib
            zgemm_(&transa,
                   &transb,
                   &nw2,
                   &nw1,
                   &nRow,
                   &alpha,
                   &wfc[ibStart][iw2_lo],
                   &gt.lgd,
                   WFC_PHASE,
                   &nw1,
                   &beta,
                   DM,
                   &nw2);

            atom2start += nw1 * nw2;
        } // ia2
    }     // ik
    return;
}

// added by zhengdy-soc, for non-collinear case
inline void cal_DM_ATOM_nc(const Grid_Technique& gt,
                           const std::complex<double> fac,
                           const Record_adj& RA,
                           const int ia1,
                           const int iw1_lo,
                           const int nw1,
                           const int gstart,
                           std::complex<double>*** wfc_k_grid,
                           std::complex<double>* WFC_PHASE,
                           std::complex<double>** DM_ATOM,
                           const ModuleBase::matrix& wg_in,
                           const K_Vectors& kv)
{

    if (GlobalV::NSPIN != 4)
    {
        ModuleBase::WARNING_QUIT("Local_Orbital_Charge", "NSPIN not match!");
    }

    const char transa = 'N';
    const char transb = 'T';
    const std::complex<double> alpha = 1;
    const std::complex<double> beta = 1;
    int ispin = 0;

    for (int is1 = 0; is1 < 2; is1++)
    {
        for (int is2 = 0; is2 < 2; is2++)
        {
            for (int ik = 0; ik < kv.get_nks(); ik++)
            {
                std::complex<double>** wfc = wfc_k_grid[ik];
                int atom2start = 0;

                for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
                {
                    std::complex<double>* DM = &DM_ATOM[ispin][atom2start];
                    const int T2 = RA.info[ia1][ia2][3];
                    const int I2 = RA.info[ia1][ia2][4];
                    Atom* atom2 = &GlobalC::ucell.atoms[T2];
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                    const int iw2_lo = gt.trace_lo[start2] / GlobalV::NPOL + gt.lgd / GlobalV::NPOL * is2;
                    const int nw2 = atom2->nw;
                    std::complex<double> exp_R = ModuleBase::libm::exp(fac
                                                                       * (kv.kvec_d[ik].x * RA.info[ia1][ia2][0]
                                                                          + kv.kvec_d[ik].y * RA.info[ia1][ia2][1]
                                                                          + kv.kvec_d[ik].z * RA.info[ia1][ia2][2]));

                    // ModuleBase::GlobalFunc::ZEROS(WFC_PHASE, GlobalV::NBANDS*nw1);
                    int ibStart = 0;
                    int nRow = 0;
                    for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
                    {
                        const double w1 = wg_in(ik, ib);
                        if (w1 > 0)
                        {
                            if (nRow == 0)
                            {
                                ibStart = ib;
                            }
                            const int iline = nRow * nw1;
                            std::complex<double> phase = exp_R * w1;
                            for (int iw1 = 0; iw1 < nw1; ++iw1)
                            {
                                WFC_PHASE[iline + iw1]
                                    = phase * conj(wfc[ib][iw1_lo + iw1 + gt.lgd / GlobalV::NPOL * is1]);
                            }
                            ++nRow;
                        }
                        else
                            break;
                    } // ib
                    zgemm_(&transa,
                           &transb,
                           &nw2,
                           &nw1,
                           &nRow,
                           &alpha,
                           &wfc[ibStart][iw2_lo],
                           &gt.lgd,
                           WFC_PHASE,
                           &nw1,
                           &beta,
                           DM,
                           &nw2);

                    atom2start += nw1 * nw2;
                } // ia2
            }     // ik
            ispin++;
        } // is2
    }     // is1
    return;
}
