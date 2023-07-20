#include "module_base/blas_connector.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/parallel_common.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "local_orbital_charge.h"
#include "module_base/libm/libm.h"

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
        this->DM_R = new double *[GlobalV::NSPIN];
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
inline void cal_DM_ATOM(const Grid_Technique &gt,
                        const std::complex<double> fac,
                        const Record_adj &RA,
                        const int ia1,
                        const int iw1_lo,
                        const int nw1,
                        const int gstart,
                        std::complex<double> ***wfc_k_grid,
                        std::complex<double> *WFC_PHASE,
                        std::complex<double> **DM_ATOM,
                        const ModuleBase::matrix &wg_in,
                        const K_Vectors& kv)
{

    const char transa = 'N';
    const char transb = 'T';
    const std::complex<double> alpha = 1;
    const std::complex<double> beta = 1;

    for (int ik = 0; ik < kv.nks; ik++)
    {
        std::complex<double> **wfc = wfc_k_grid[ik];
        const int ispin = kv.isk[ik];
        int atom2start = 0;

        for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
        {
            std::complex<double> *DM = &DM_ATOM[ispin][atom2start];
            const int T2 = RA.info[ia1][ia2][3];
            const int I2 = RA.info[ia1][ia2][4];
            Atom *atom2 = &GlobalC::ucell.atoms[T2];
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
    } // ik
    return;
}

// added by zhengdy-soc, for non-collinear case
inline void cal_DM_ATOM_nc(const Grid_Technique &gt,
                           const std::complex<double> fac,
                           const Record_adj &RA,
                           const int ia1,
                           const int iw1_lo,
                           const int nw1,
                           const int gstart,
                           std::complex<double> ***wfc_k_grid,
                           std::complex<double> *WFC_PHASE,
                           std::complex<double> **DM_ATOM,
                           const ModuleBase::matrix &wg_in,
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
            for (int ik = 0; ik < kv.nks; ik++)
            {
                std::complex<double> **wfc = wfc_k_grid[ik];
                int atom2start = 0;

                for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
                {
                    std::complex<double> *DM = &DM_ATOM[ispin][atom2start];
                    const int T2 = RA.info[ia1][ia2][3];
                    const int I2 = RA.info[ia1][ia2][4];
                    Atom *atom2 = &GlobalC::ucell.atoms[T2];
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
            } // ik
            ispin++;
        } // is2
    } // is1
    return;
}

void Local_Orbital_Charge::cal_dk_k(const Grid_Technique &gt, const ModuleBase::matrix &wg_in, const K_Vectors& kv)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "cal_dk_k");
    ModuleBase::timer::tick("LCAO_Charge", "cal_dk_k");
    // int nnrg = 0;

    Record_adj RA;
    RA.for_grid(gt);

#ifdef __MKL
		const int mkl_threads = mkl_get_max_threads();
		mkl_set_num_threads(1);
#endif

#ifdef _OPENMP
#pragma omp parallel
{
#endif
    ModuleBase::Vector3<double> tau1;
    ModuleBase::Vector3<double> dtau;
    std::complex<double> fac = ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT;

    std::complex<double> *WFC_PHASE = new std::complex<double>[GlobalV::NLOCAL * GlobalC::ucell.nwmax];

    int DM_ATOM_SIZE = 1;
    std::complex<double> **DM_ATOM = new std::complex<double> *[GlobalV::NSPIN];

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        DM_ATOM[is] = new std::complex<double>[DM_ATOM_SIZE];
        ModuleBase::GlobalFunc::ZEROS(DM_ATOM[is], DM_ATOM_SIZE);
    }
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
    for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
	{
		const int T1 = GlobalC::ucell.iat2it[iat];
		Atom* atom1 = &GlobalC::ucell.atoms[T1];
		const int I1 = GlobalC::ucell.iat2ia[iat];
		{
			const int ca = RA.iat2ca[iat];
            if (gt.in_this_processor[iat])
            {
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const int gstart = gt.nlocstartg[iat];
                const int ng = gt.nlocdimg[iat];
                const int iw1_lo = gt.trace_lo[start1] / GlobalV::NPOL;
                const int nw1 = atom1->nw;

                if (DM_ATOM_SIZE < ng)
                {
                    DM_ATOM_SIZE = ng;
                    for (int is = 0; is < GlobalV::NSPIN; ++is)
                    {
                        delete[] DM_ATOM[is];
                    }
                    for (int is = 0; is < GlobalV::NSPIN; ++is)
                    {
                        DM_ATOM[is] = new std::complex<double>[DM_ATOM_SIZE];
                    }
                }
                for (int is = 0; is < GlobalV::NSPIN; ++is)
                {
                    ModuleBase::GlobalFunc::ZEROS(DM_ATOM[is], ng);
                }
                ModuleBase::GlobalFunc::ZEROS(WFC_PHASE, GlobalV::NBANDS * nw1);
                if (GlobalV::NSPIN != 4)
                {
                    cal_DM_ATOM(gt,
                                fac,
                                RA,
                                ca,
                                iw1_lo,
                                nw1,
                                gstart,
                                this->LOWF->wfc_k_grid,
                                WFC_PHASE,
                                DM_ATOM,
                                wg_in,
                                kv);
                }
                else
                {
                    cal_DM_ATOM_nc(gt,
                                   fac,
                                   RA,
                                   ca,
                                   iw1_lo,
                                   nw1,
                                   gstart,
                                   this->LOWF->wfc_k_grid,
                                   WFC_PHASE,
                                   DM_ATOM,
                                   wg_in,
                                   kv);
                }

                if (GlobalV::NSPIN != 4)
                {
                    for (int is = 0; is < GlobalV::NSPIN; ++is)
                    {
                        for (int iv = 0; iv < ng; ++iv)
                        {
                            this->DM_R[is][gstart + iv] = DM_ATOM[is][iv].real();
                        }
                    }
                }
                else
                { // zhengdy-soc
                    for (int iv = 0; iv < ng; ++iv)
                    {
                        // note: storage nondiagonal term as Re[] and Im[] respectly;
                        this->DM_R[0][gstart + iv] = DM_ATOM[0][iv].real() + DM_ATOM[3][iv].real();
                        if (GlobalV::NONCOLIN)
                        { // GlobalV::DOMAG
                            this->DM_R[1][gstart + iv] = DM_ATOM[1][iv].real() + DM_ATOM[2][iv].real();
                            this->DM_R[2][gstart + iv] = DM_ATOM[1][iv].imag() - DM_ATOM[2][iv].imag();
                            this->DM_R[3][gstart + iv] = DM_ATOM[0][iv].real() - DM_ATOM[3][iv].real();
                        }
                        else if (!GlobalV::NONCOLIN) // GlobalV::DOMAG_Z
                        {
                            this->DM_R[1][gstart + iv] = 0.0;
                            this->DM_R[2][gstart + iv] = 0.0;
                            this->DM_R[3][gstart + iv] = DM_ATOM[0][iv].real() - DM_ATOM[3][iv].real();
                        }
                        else // soc with no mag
                        {
                            this->DM_R[1][gstart + iv] = 0.0;
                            this->DM_R[2][gstart + iv] = 0.0;
                            this->DM_R[3][gstart + iv] = 0.0;
                        }
                    }
                }
            } // if gt.in_this_processor
        } // I1
    } // T1

    //------------
    // for test
    //------------
    /*  std::cout << std::setprecision(3);
        for(int i=0; i<nnrg_now; i++)

        for(int ik=0; ik<kv.nkstot; ++ik)
        {
            for(int ib=0; ib<GlobalV::NBANDS; ++ib)
            {
                std::cout << " ik=" << ik << " ib=" << ib << " occ=" << GlobalC::wf.wg(ik,ib) << " e=" <<
       GlobalC::wf.ekb[ik][ib] << std::endl;
            }
        }

        for(int i=0; i<10; i++)
        {
            if(DM_R[0][i]>1.0e-8)
            {
                std::cout << " i=" << i << " DM_R=" << DM_R[0][i] << std::endl;
            }
        }
    */
    for (int i = 0; i < GlobalV::NSPIN; ++i)
    {
        delete[] DM_ATOM[i];
    }
    delete[] DM_ATOM;
    delete[] WFC_PHASE;
#ifdef _OPENMP
}
#endif

#ifdef __MKL
    mkl_set_num_threads(mkl_threads);
#endif

    RA.delete_grid(); // xiaohui add 2015-02-04

    ModuleBase::timer::tick("LCAO_Charge", "cal_dk_k");
    return;
}
