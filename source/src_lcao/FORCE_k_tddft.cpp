#include "FORCE_k_tddft.h"

#include "../module_base/memory.h"
#include "../module_base/timer.h"
#include "../src_parallel/parallel_reduce.h"
#include "../src_pw/global.h"
#include "ELEC_evolve.h"
#include "module_elecstate/cal_dm.h"

#include <map>
#include <unordered_map>

#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif

Force_LCAO_k_tddft::Force_LCAO_k_tddft()
{
}

Force_LCAO_k_tddft::~Force_LCAO_k_tddft()
{
}

// be called in Force_LCAO::start_force_calculation
void Force_LCAO_k_tddft::ftable_k(const bool isforce,
                                  const bool isstress,
                                  Record_adj& ra,
                                  const psi::Psi<std::complex<double>>* psi,
                                  Local_Orbital_Charge& loc,
                                  ModuleBase::matrix& foverlap,
                                  ModuleBase::matrix& ftddft,
                                  ModuleBase::matrix& ftvnl_dphi,
                                  ModuleBase::matrix& fvnl_dbeta,
                                  ModuleBase::matrix& fvl_dphi,
                                  ModuleBase::matrix& soverlap,
                                  ModuleBase::matrix& stvnl_dphi,
                                  ModuleBase::matrix& svnl_dbeta,
#ifdef __DEEPKS
                                  ModuleBase::matrix& svl_dphi,
                                  ModuleBase::matrix& svnl_dalpha,
#else
                                  ModuleBase::matrix& svl_dphi,
#endif
                                  LCAO_Hamilt& uhm,
                                  ModuleBase::Vector3<double>* vel)
{
    ModuleBase::TITLE("Force_LCAO_k_tddft", "ftable_k");
    ModuleBase::timer::tick("Force_LCAO_k_tddft", "ftable_k");

    this->UHM = &uhm;

    const Parallel_Orbitals* pv = loc.ParaV;
    this->allocate_k(*pv);

    // calculate the energy density matrix
    // and the force related to overlap matrix and energy density matrix.
    // Force_LCAO_k::cal_foverlap_k(isforce, isstress, ra, psi, loc, foverlap, soverlap);
    this->cal_foverlap_k(isforce, isstress, ra, psi, loc, foverlap, soverlap);

    // calculate the density matrix
    double** dm2d = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        dm2d[is] = new double[pv->nnr];
        ModuleBase::GlobalFunc::ZEROS(dm2d[is], pv->nnr);
    }
    ModuleBase::Memory::record("Force_LCAO_k", "dm2d", GlobalV::NSPIN * pv->nnr, "double");

    loc.cal_dm_R(loc.dm_k, ra, dm2d);

    this->cal_ftvnl_dphi_k(dm2d, isforce, isstress, ra, ftvnl_dphi, stvnl_dphi);

    //this->cal_ftddft_k(isforce, ra, psi, loc, ftddft, vel);

    // ---------------------------------------
    // doing on the real space grid.
    // ---------------------------------------
    this->cal_fvl_dphi_k(isforce, isstress, fvl_dphi, svl_dphi, loc.DM_R);

    this->calFvnlDbeta(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta, GlobalV::vnl_method);

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        delete[] dm2d[is];
    }
    delete[] dm2d;

    //----------------------------------------------------------------
    // reduce the force according to 2D distribution of H & S matrix.
    //----------------------------------------------------------------
    if (isforce)
    {
        Parallel_Reduce::reduce_double_pool(foverlap.c, foverlap.nr * foverlap.nc);
        Parallel_Reduce::reduce_double_pool(ftvnl_dphi.c, ftvnl_dphi.nr * ftvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool(fvnl_dbeta.c, fvnl_dbeta.nr * fvnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool(fvl_dphi.c, fvl_dphi.nr * fvl_dphi.nc);
    }

    this->finish_k();

    ModuleBase::timer::tick("Force_LCAO_k", "ftable_k");
    return;
}
///*
void Force_LCAO_k_tddft::allocate_k(const Parallel_Orbitals& pv)
{
    ModuleBase::TITLE("Force_LCAO_k_tddft", "allocate_k");
    ModuleBase::timer::tick("Force_LCAO_k_tddft", "allocate_k");

    this->ParaV = &pv;
    const int nnr = pv.nnr;
    //--------------------------------
    // (1) allocate for dSx dSy & dSz
    //--------------------------------
    this->UHM->LM->DSloc_Rx = new double[nnr];
    this->UHM->LM->DSloc_Ry = new double[nnr];
    this->UHM->LM->DSloc_Rz = new double[nnr];
    ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DSloc_Rx, nnr);
    ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DSloc_Ry, nnr);
    ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DSloc_Rz, nnr);
    ModuleBase::Memory::record("force_lo", "dS", nnr * 3, "double");

    if (GlobalV::CAL_STRESS)
    {
        this->UHM->LM->DH_r = new double[3 * nnr];
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DH_r, 3 * nnr);
        this->UHM->LM->stvnl11 = new double[nnr];
        this->UHM->LM->stvnl12 = new double[nnr];
        this->UHM->LM->stvnl13 = new double[nnr];
        this->UHM->LM->stvnl22 = new double[nnr];
        this->UHM->LM->stvnl23 = new double[nnr];
        this->UHM->LM->stvnl33 = new double[nnr];
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl11, nnr);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl12, nnr);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl13, nnr);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl22, nnr);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl23, nnr);
        ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->stvnl33, nnr);
        ModuleBase::Memory::record("stress_lo", "dSR", nnr * 6, "double");
    }

    //-----------------------------
    // calculate dS = <phi | dphi>
    //-----------------------------
    bool cal_deri = true;
    this->UHM->genH.build_ST_new('S', cal_deri, GlobalC::ucell, this->UHM->genH.LM->SlocR.data());

    //-----------------------------------------
    // (2) allocate for <phi | T + Vnl | dphi>
    //-----------------------------------------
    this->UHM->LM->DHloc_fixedR_x = new double[nnr];
    this->UHM->LM->DHloc_fixedR_y = new double[nnr];
    this->UHM->LM->DHloc_fixedR_z = new double[nnr];
    ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixedR_x, nnr);
    ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixedR_y, nnr);
    ModuleBase::GlobalFunc::ZEROS(this->UHM->LM->DHloc_fixedR_z, nnr);
    ModuleBase::Memory::record("force_lo", "dTVNL", nnr * 3, "double");

    // calculate dT=<phi|kin|dphi> in LCAO
    // calculate T + VNL(P1) in LCAO basis
    this->UHM->genH.build_ST_new('T', cal_deri, GlobalC::ucell, this->UHM->genH.LM->Hloc_fixedR.data());
    // test(this->UHM->LM->DHloc_fixedR_x,"this->UHM->LM->DHloc_fixedR_x T part");

    // calculate dVnl=<phi|dVnl|dphi> in LCAO
    this->NonlocalDphi(GlobalV::NSPIN, GlobalV::vnl_method, cal_deri, this->UHM->genH);
    // test(this->UHM->LM->DHloc_fixedR_x,"this->UHM->LM->DHloc_fixedR_x Vnl part");

    ModuleBase::timer::tick("Force_LCAO_k", "allocate_k");
    return;
} //*/

void Force_LCAO_k_tddft::finish_k(void)
{
    delete[] this->UHM->LM->DSloc_Rx;
    delete[] this->UHM->LM->DSloc_Ry;
    delete[] this->UHM->LM->DSloc_Rz;
    /*for (int i = 0; i < this->ParaV->nnr; i++)
    {
        delete[] this->UHM->LM->DDSloc[i];
    }
    delete[] this->UHM->LM->DDSloc;*/
    delete[] this->UHM->LM->DHloc_fixedR_x;
    delete[] this->UHM->LM->DHloc_fixedR_y;
    delete[] this->UHM->LM->DHloc_fixedR_z;
    if (GlobalV::CAL_STRESS)
    {
        delete[] this->UHM->LM->DH_r;
        delete[] this->UHM->LM->stvnl11;
        delete[] this->UHM->LM->stvnl12;
        delete[] this->UHM->LM->stvnl13;
        delete[] this->UHM->LM->stvnl22;
        delete[] this->UHM->LM->stvnl23;
        delete[] this->UHM->LM->stvnl33;
    }
    return;
}

#include "record_adj.h"
void Force_LCAO_k_tddft::cal_ftddft_k(const bool isforce,
                                      Record_adj& ra,
                                      const psi::Psi<std::complex<double>>* psi,
                                      Local_Orbital_Charge& loc,
                                      ModuleBase::matrix& ftddft,
                                      ModuleBase::Vector3<double>* vel)
{
    ModuleBase::TITLE("Force_LCAO_k_tddft", "cal_ftddft_k");
    ModuleBase::timer::tick("Force_LCAO_k_tddft", "cal_ftddft_k");

    const Parallel_Orbitals* pv = this->ParaV;

    complex<double> C1x, C2x, C1y, C2y, C1z, C2z;

    complex<double> imag = {0, 1.0};
    int nrow = this->ParaV->nrow;
    int ncol = this->ParaV->ncol;
    int ncol_bands = this->ParaV->ncol_bands;
    int nloc = this->ParaV->nloc;

    complex<double>* Sinv = new complex<double>[nloc];
    const int inc = 1;
    ModuleBase::GlobalFunc::ZEROS(Sinv, nloc);
    zcopy_(&this->ParaV->nloc, this->UHM->LM->Sloc2.data(), &inc, Sinv, &inc);

    int* ipiv = new int[nloc];
    int info;
    const int one_int = 1;
    pzgetrf_(&GlobalV::NLOCAL, &GlobalV::NLOCAL, Sinv, &one_int, &one_int, this->ParaV->desc, ipiv, &info);

    int LWORK = -1, liWORK = -1;
    std::vector<std::complex<double>> WORK(1, 0);
    std::vector<int> iWORK(1, 0);

    GlobalV::ofs_running << "test1" << endl;
    pzgetri_(&GlobalV::NLOCAL,
             Sinv,
             &one_int,
             &one_int,
             this->ParaV->desc,
             ipiv,
             WORK.data(),
             &LWORK,
             iWORK.data(),
             &liWORK,
             &info);

    LWORK = WORK[0].real();
    WORK.resize(LWORK, 0);
    liWORK = iWORK[0];
    iWORK.resize(liWORK, 0);

    pzgetri_(&GlobalV::NLOCAL,
             Sinv,
             &one_int,
             &one_int,
             this->ParaV->desc,
             ipiv,
             WORK.data(),
             &LWORK,
             iWORK.data(),
             &liWORK,
             &info);

    for (int i = 0; i < GlobalC::ucell.nat; i++)
    {
        ftddft(i, 0) = 0.0;
        ftddft(i, 1) = 0.0;
        ftddft(i, 2) = 0.0;
    }

    // GlobalV::ofs_running << endl;
    // GlobalV::ofs_running << "print initial ftddft: " << endl;
    // for (int i = 0; i < GlobalC::ucell.nat; i++)
    // {
    //     GlobalV::ofs_running << i << " " << ftddft(i, 0) << " " << ftddft(i, 1) << " " << ftddft(i, 2) << endl;
    // }
    // GlobalV::ofs_running << endl;

    for (int ik = 0; ik < GlobalC::kv.nks; ik++)
    {
        psi->fix_k(ik);
        int irr1 = 0;
        int iat1 = 0;
        for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
        {
            Atom* atom1 = &GlobalC::ucell.atoms[T1];
            for (int I1 = 0; I1 < atom1->na; ++I1)
            {
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                int iat2 = 0;
                for (int cb = 0; cb < ra.na_each[iat1]; ++cb)
                {
                    const int T2 = ra.info[iat1][cb][3];
                    const int I2 = ra.info[iat1][cb][4];
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

                    Atom* atom2 = &GlobalC::ucell.atoms[T2];

                    for (int jj = 0; jj < atom1->nw; jj++)
                    {
                        const int iw1_all = start1 + jj;

                        const int mu = pv->trace_loc_row[iw1_all];
                        if (mu < 0)
                            continue;

                        for (int kk = 0; kk < atom2->nw; kk++)
                        {
                            const int iw2_all = start2 + kk;

                            const int nu = pv->trace_loc_col[iw2_all];
                            if (nu < 0)
                                continue;

                            int irr2 = 0;
                            int iat3 = 0;
                            for (int T3 = 0; T3 < GlobalC::ucell.ntype; ++T3)
                            {
                                Atom* atom3 = &GlobalC::ucell.atoms[T3];
                                for (int I3 = 0; I3 < atom3->na; ++I3)
                                {
                                    const int start3 = GlobalC::ucell.itiaiw2iwt(T3, I3, 0);
                                    int iat4 = 0;
                                    for (int cbb = 0; cbb < ra.na_each[iat3]; ++cbb)
                                    {
                                        const int T4 = ra.info[iat3][cbb][3];
                                        const int I4 = ra.info[iat3][cbb][4];
                                        const int start4 = GlobalC::ucell.itiaiw2iwt(T4, I4, 0);

                                        Atom* atom4 = &GlobalC::ucell.atoms[T4];

                                        for (int ll = 0; ll < atom3->nw; ll++)
                                        {
                                            const int iw3_all = start3 + ll;

                                            const int eta = pv->trace_loc_row[iw3_all];
                                            if (eta < 0)
                                                continue;

                                            for (int mm = 0; mm < atom4->nw; mm++)
                                            {
                                                const int iw4_all = start4 + mm;

                                                const int theta = pv->trace_loc_col[iw4_all];
                                                if (theta < 0)
                                                    continue;
                                                for (int iband = 0; iband < ncol_bands; ++iband)
                                                {
                                                    complex<double> wfc_part, B_part1_x, B_part1_y, B_part1_z,
                                                        B_part2_x, B_part2_y, B_part2_z;
                                                    wfc_part = conj(psi[0].get_pointer()[iband * ncol + mu])
                                                               * psi[0].get_pointer()[iband * ncol + eta];
                                                    // GlobalV::ofs_running << "wfc_part= "<<wfc_part <<endl;
                                                    B_part1_x = vel[iat1][0] * this->UHM->LM->DSloc_Rx[irr1]
                                                                * Sinv[nu * ncol + theta]
                                                                * this->UHM->LM->DSloc_Rx[irr2] * wfc_part;
                                                    B_part1_y = vel[iat1][1] * this->UHM->LM->DSloc_Ry[irr1]
                                                                * Sinv[nu * ncol + theta]
                                                                * this->UHM->LM->DSloc_Ry[irr2] * wfc_part;
                                                    B_part1_z = vel[iat1][2] * this->UHM->LM->DSloc_Rz[irr1]
                                                                * Sinv[nu * ncol + theta]
                                                                * this->UHM->LM->DSloc_Rz[irr2] * wfc_part;
                                                    B_part2_x = -vel[iat3][0] * this->UHM->LM->DSloc_Rx[irr1]
                                                                * Sinv[nu * ncol + theta]
                                                                * this->UHM->LM->DSloc_Rx[irr2] * wfc_part;
                                                    B_part2_y = -vel[iat3][1] * this->UHM->LM->DSloc_Ry[irr1]
                                                                * Sinv[nu * ncol + theta]
                                                                * this->UHM->LM->DSloc_Ry[irr2] * wfc_part;
                                                    B_part2_z = -vel[iat3][2] * this->UHM->LM->DSloc_Rz[irr1]
                                                                * Sinv[nu * ncol + theta]
                                                                * this->UHM->LM->DSloc_Rz[irr2] * wfc_part;
                                                    //GlobalV::ofs_running << "vel[iat1]: "<<vel[iat1][0]<<" "<<vel[iat1][1]<<" "<<vel[iat1][2]<<endl;
                                                    //GlobalV::ofs_running << "vel[iat3]: "<<vel[iat3][0]<<" "<<vel[iat3][1]<<" "<<vel[iat3][2]<<endl;
                                                    //GlobalV::ofs_running << "DS[irr1]="<<this->UHM->LM->DSloc_Ry[irr1]<<" DS[irr2]="<<this->UHM->LM->DSloc_Ry[irr2] <<endl;
                                                    // GlobalV::ofs_running<<"iat : "<<iat1<<" "<<iat3<<endl;
                                                    // cout<<"size of ftddft"<<sizeof(ftddft)<<endl;
                                                    /*ftddft(iat1, 0) = 0.0;
                                                    ftddft(iat1, 1) = 0.0;
                                                    ftddft(iat1, 2) = 0.0;*/

                                                    //GlobalV::ofs_running << "iat1 = " << iat1 << " iat3= " << iat3
                                                    //                     << endl;

                                                    // GlobalV::ofs_running << "ftddft iat1 before " << ftddft(iat1, 0)
                                                    //                      << " " << ftddft(iat1, 1) << " "
                                                    //                      << ftddft(iat1, 2) << endl;
                                                    // GlobalV::ofs_running << "ftddft iat3 before " << ftddft(iat3, 0)
                                                    //                      << " " << ftddft(iat3, 1) << " "
                                                    //                      << ftddft(iat3, 2) << endl;

                                                    ftddft(iat3, 0) += B_part1_x.imag();
                                                    ftddft(iat3, 1) += B_part1_y.imag();
                                                    ftddft(iat3, 2) += B_part1_z.imag();
                                                    ftddft(iat1, 0) += B_part2_x.imag();
                                                    ftddft(iat1, 1) += B_part2_y.imag();
                                                    ftddft(iat1, 2) += B_part2_z.imag();

                                                    // GlobalV::ofs_running << "1x :" << B_part1_x.imag() << endl;
                                                    // GlobalV::ofs_running << "1y :" << B_part1_y.imag() << endl;
                                                    // GlobalV::ofs_running << "1z :" << B_part1_z.imag() << endl;
                                                    // GlobalV::ofs_running << "2x :" << B_part2_x.imag() << endl;
                                                    // GlobalV::ofs_running << "2y :" << B_part2_y.imag() << endl;
                                                    // GlobalV::ofs_running << "2z :" << B_part2_z.imag() << endl;

                                                    // GlobalV::ofs_running << "ftddft iat1 after " << ftddft(iat1, 0)
                                                    //                      << " " << ftddft(iat1, 1) << " "
                                                    //                      << ftddft(iat1, 2) << endl;
                                                    // GlobalV::ofs_running << "ftddft iat3 after " << ftddft(iat3, 0)
                                                    //                      << " " << ftddft(iat3, 1) << " "
                                                    //                      << ftddft(iat3, 2) << endl;

                                                    // GlobalV::ofs_running << "iat1=" << iat1 << " iat3=" << iat3 <<
                                                    // endl;
                                                    // double tmp;
                                                    // tmp = 0.01;
                                                    // if (abs(B_part1_x.imag()) > tmp)
                                                    // {
                                                    //     GlobalV::ofs_running << "1x :" << B_part1_x.imag() << endl;
                                                    // }
                                                    // if (abs(B_part1_y.imag()) > tmp)
                                                    // {
                                                    //     GlobalV::ofs_running << "1y :" << B_part1_y.imag() << endl;
                                                    // }
                                                    // if (abs(B_part1_z.imag()) > tmp)
                                                    // {
                                                    //     GlobalV::ofs_running << "1z :" << B_part1_z.imag() << endl;
                                                    // }
                                                    // if (abs(B_part2_x.imag()) > tmp)
                                                    // {
                                                    //     GlobalV::ofs_running << "2x :" << B_part2_x.imag() << endl;
                                                    // }
                                                    // if (abs(B_part2_y.imag()) > tmp)
                                                    // {
                                                    //     GlobalV::ofs_running << "2y :" << B_part2_y.imag() << endl;
                                                    // }
                                                    // if (abs(B_part2_z.imag()) > tmp)
                                                    // {
                                                    //     GlobalV::ofs_running << "2z :" << B_part2_z.imag() << endl;
                                                    // }
                                                }
                                                ++irr2;
                                            }
                                        }
                                        ++iat4;
                                    }
                                    ++iat3;
                                }
                            }

                            ++irr1;
                        } // end kk
                    } // end jj
                    ++iat2;
                } // end cb
                ++iat1;
            }
        }
    }

    GlobalV::ofs_running << endl;
    GlobalV::ofs_running << "print ftddft: " << endl;
    for (int i = 0; i < GlobalC::ucell.nat; i++)
    {
        GlobalV::ofs_running << i << " " << ftddft(i, 0) << " " << ftddft(i, 1) << " " << ftddft(i, 2) << endl;
    }
    GlobalV::ofs_running << endl;

    delete[] Sinv;
    /*
    if (irr1 != pv->nnr || irr2 != pv->nnr)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong irr", irr);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong LNNR.nnr", pv->nnr);
        ModuleBase::WARNING_QUIT("Force_LCAO_k_tddft::cal_ftddft_k", "irr!=LNNR.nnr");
    }*/

    ModuleBase::timer::tick("Force_LCAO_k_tddft", "cal_ftddft_k");
    return;
}

/*
#include "record_adj.h"
void Force_LCAO_k_tddft::cal_foverlap_k(const bool isforce,
                                  const bool isstress,
                                  Record_adj& ra,
                                  const psi::Psi<std::complex<double>>* psi,
                                  Local_Orbital_Charge& loc,
                                  ModuleBase::matrix& foverlap,
                                  ModuleBase::matrix& soverlap)
{
    ModuleBase::TITLE("Force_LCAO_k_tddft", "cal_foverlap_k");
    ModuleBase::timer::tick("Force_LCAO_k_tddft", "cal_foverlap_k");

    const Parallel_Orbitals* pv = this->ParaV;
    //--------------------------------------------
    // (1) allocate energy density matrix (nnr)
    //--------------------------------------------
    double** edm2d = new double*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        edm2d[is] = new double[pv->nnr];
        ModuleBase::GlobalFunc::ZEROS(edm2d[is], pv->nnr);
    }

    //--------------------------------------------
    // calculate the energy density matrix here.
    //--------------------------------------------
    ModuleBase::timer::tick("Force_LCAO_k", "cal_edm_2d");

    ModuleBase::matrix wgEkb;
    wgEkb.create(GlobalC::kv.nks, GlobalV::NBANDS);

    for (int ik = 0; ik < GlobalC::kv.nks; ik++)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            wgEkb(ik, ib) = GlobalC::wf.wg(ik, ib) * GlobalC::wf.ekb[ik][ib];
        }
    }
    std::vector<ModuleBase::ComplexMatrix> edm_k(GlobalC::kv.nks);

    // use the original formula (Hamiltonian matrix) to calculate energy density matrix
    if (loc.edm_k_tddft.size())
    {
        GlobalV::ofs_running<<"td edm"<<endl;
        for (int ik = 0; ik < GlobalC::kv.nks; ++ik)
        {
            edm_k[ik] = loc.edm_k_tddft[ik];
        }
    }
    else
    {
        elecstate::cal_dm(loc.ParaV, wgEkb, psi[0], edm_k);
    }

    GlobalV::ofs_running<<"before cal_dm_R trace_loc_row="<<pv->trace_loc_row<<endl;

    loc.cal_dm_R(edm_k, ra, edm2d);
    ModuleBase::timer::tick("Force_LCAO_k", "cal_edm_2d");

    GlobalV::ofs_running<<"after cal_dm_R trace_loc_row="<<sizeof(pv->trace_loc_row)<<endl;

    //--------------------------------------------
    // summation \sum_{i,j} E(i,j)*dS(i,j)
    // BEGIN CALCULATION OF FORCE OF EACH ATOM
    //--------------------------------------------
    ModuleBase::Vector3<double> tau1, dtau, tau2;

    GlobalV::ofs_running<<"foverlap test0"<<endl;
    int irr = 0;
    int iat = 0;
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            for (int cb = 0; cb < ra.na_each[iat]; ++cb)
            {
                const int T2 = ra.info[iat][cb][3];
                const int I2 = ra.info[iat][cb][4];
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);

                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                for (int jj = 0; jj < atom1->nw; jj++)
                {
                    const int iw1_all = start1 + jj;
                    GlobalV::ofs_running<<"foverlap test01"<<endl;
                    GlobalV::ofs_running<<"size="<<sizeof(pv->trace_loc_row)<<endl;
                    GlobalV::ofs_running<<"iw1_all="<<iw1_all<<endl;
                    GlobalV::ofs_running<<"output 1"<<pv->trace_loc_row[1]<<endl;
                    GlobalV::ofs_running<<"output 0"<<pv->trace_loc_row[0]<<endl;
                    // HPSEPS
                    const int mu = pv->trace_loc_row[iw1_all];
                    GlobalV::ofs_running<<"trace_loc_row="<<pv->trace_loc_row<<endl;
                    //GlobalV::ofs_running<<"trace_loc_row[iw1_all]="<<pv->trace_loc_row[iw1_all]<<endl;
                    GlobalV::ofs_running<<"foverlap test03"<<endl;
                    if (mu < 0)
                        continue;
                    GlobalV::ofs_running<<"foverlap test04"<<endl;
                    for (int kk = 0; kk < atom2->nw; kk++)
                    {
                        GlobalV::ofs_running<<"foverlap test05"<<endl;
                        const int iw2_all = start2 + kk;
                        GlobalV::ofs_running<<"foverlap test02"<<endl;
                        // HPSEPS
                        const int nu = pv->trace_loc_col[iw2_all];
                        if (nu < 0)
                            continue;
                        //==============================================================
                        // here we use 'minus', but in GlobalV::GAMMA_ONLY_LOCAL we use 'plus',
                        // both are correct because the 'DSloc_Rx' is used in 'row' (-),
                        // however, the 'DSloc_x' in GAMMA is used in 'col' (+),
                        // mohan update 2011-06-16
                        //==============================================================
                        for (int is = 0; is < GlobalV::NSPIN; ++is)
                        {
                            GlobalV::ofs_running<<"foverlap test1"<<endl;
                            double edm2d2 = 2.0 * edm2d[is][irr];
                            GlobalV::ofs_running<<"foverlap test2"<<endl;
                            if (isforce)
                            {
                                foverlap(iat, 0) -= edm2d2 * this->UHM->LM->DSloc_Rx[irr];
                                foverlap(iat, 1) -= edm2d2 * this->UHM->LM->DSloc_Ry[irr];
                                foverlap(iat, 2) -= edm2d2 * this->UHM->LM->DSloc_Rz[irr];
                            }
                            if (isstress)
                            {
                                for (int ipol = 0; ipol < 3; ipol++)
                                {
                                    soverlap(0, ipol) += edm2d[is][irr] * this->UHM->LM->DSloc_Rx[irr]
                                                         * this->UHM->LM->DH_r[irr * 3 + ipol];
                                    if (ipol < 1)
                                        continue;
                                    soverlap(1, ipol) += edm2d[is][irr] * this->UHM->LM->DSloc_Ry[irr]
                                                         * this->UHM->LM->DH_r[irr * 3 + ipol];
                                    if (ipol < 2)
                                        continue;
                                    soverlap(2, ipol) += edm2d[is][irr] * this->UHM->LM->DSloc_Rz[irr]
                                                         * this->UHM->LM->DH_r[irr * 3 + ipol];
                                }
                            }
                        }
                        ++irr;
                    } // end kk
                } // end jj
            } // end cb
            ++iat;
        }
    }

    GlobalV::ofs_running<<"foverlap"<<endl;

    if (isstress)
    {
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, soverlap);
    }

    if (irr != pv->nnr)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong irr", irr);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "wrong LNNR.nnr", pv->nnr);
        ModuleBase::WARNING_QUIT("Force_LCAO_k::cal_foverlap_k", "irr!=LNNR.nnr");
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        delete[] edm2d[is];
    }
    delete[] edm2d;

    ModuleBase::timer::tick("Force_LCAO_k_tddft", "cal_foverlap_k");
    return;
}

void Force_LCAO_k_tddft::cal_ftvnl_dphi_k(double** dm2d,
                                    const bool isforce,
                                    const bool isstress,
                                    Record_adj& ra,
                                    ModuleBase::matrix& ftvnl_dphi,
                                    ModuleBase::matrix& stvnl_dphi)
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_ftvnl_dphi");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_ftvnl_dphi");

    const Parallel_Orbitals* pv = this->ParaV;
    // get the adjacent atom's information.

    //	GlobalV::ofs_running << " calculate the ftvnl_dphi_k force" << std::endl;

    int irr = 0;
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            const int iat = GlobalC::ucell.itia2iat(T1, I1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            for (int cb = 0; cb < ra.na_each[iat]; ++cb)
            {
                const int T2 = ra.info[iat][cb][3];
                const int I2 = ra.info[iat][cb][4];
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                Atom* atom2 = &GlobalC::ucell.atoms[T2];

                for (int jj = 0; jj < atom1->nw; ++jj)
                {
                    const int iw1_all = start1 + jj;
                    const int mu = pv->trace_loc_row[iw1_all];
                    if (mu < 0)
                        continue;
                    for (int kk = 0; kk < atom2->nw; ++kk)
                    {
                        const int iw2_all = start2 + kk;
                        const int nu = pv->trace_loc_col[iw2_all];
                        if (nu < 0)
                            continue;
                        //==============================================================
                        // here we use 'minus', but in GlobalV::GAMMA_ONLY_LOCAL we use 'plus',
                        // both are correct because the 'DSloc_Rx' is used in 'row' (-),
                        // however, the 'DSloc_x' is used in 'col' (+),
                        // mohan update 2011-06-16
                        //==============================================================
                        for (int is = 0; is < GlobalV::NSPIN; ++is)
                        {
                            double dm2d2 = 2.0 * dm2d[is][irr];
                            if (isforce)
                            {
                                ftvnl_dphi(iat, 0) += dm2d2 * this->UHM->LM->DHloc_fixedR_x[irr];
                                ftvnl_dphi(iat, 1) += dm2d2 * this->UHM->LM->DHloc_fixedR_y[irr];
                                ftvnl_dphi(iat, 2) += dm2d2 * this->UHM->LM->DHloc_fixedR_z[irr];
                            }
                            if (isstress)
                            {
                                stvnl_dphi(0, 0) -= dm2d[is][irr] * this->UHM->LM->stvnl11[irr];
                                stvnl_dphi(0, 1) -= dm2d[is][irr] * this->UHM->LM->stvnl12[irr];
                                stvnl_dphi(0, 2) -= dm2d[is][irr] * this->UHM->LM->stvnl13[irr];
                                stvnl_dphi(1, 1) -= dm2d[is][irr] * this->UHM->LM->stvnl22[irr];
                                stvnl_dphi(1, 2) -= dm2d[is][irr] * this->UHM->LM->stvnl23[irr];
                                stvnl_dphi(2, 2) -= dm2d[is][irr] * this->UHM->LM->stvnl33[irr];
                            }
                        }
                        ++irr;
                    } // end kk
                } // end jj
            } // end cb
        }
    }
    assert(irr == pv->nnr);

    //	test(this->UHM->LM->DSloc_Rx);
    //	test(dm2d[0],"dm2d");

    if (isstress)
    {
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, stvnl_dphi);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "cal_ftvnl_dphi");
    return;
}

// calculate the force due to < phi | Vlocal | dphi >
void Force_LCAO_k_tddft::cal_fvl_dphi_k(const bool isforce,
                                  const bool isstress,
                                  ModuleBase::matrix& fvl_dphi,
                                  ModuleBase::matrix& svl_dphi,
                                  double** DM_R)
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_fvl_dphi_k");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");

    if (!isforce && !isstress)
    {
        ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");
        return;
    }
    assert(this->UHM->LM->DHloc_fixedR_x != NULL);
    assert(this->UHM->LM->DHloc_fixedR_y != NULL);
    assert(this->UHM->LM->DHloc_fixedR_z != NULL);

    int istep = 1;

    // if Vna potential is not used.
    GlobalC::pot.init_pot(istep, GlobalC::sf.strucFac);

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;

        for (int ir = 0; ir < GlobalC::rhopw->nrxx; ir++)
        {
            GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        }

        //--------------------------------
        // Grid integration here.
        //--------------------------------
        // fvl_dphi can not be set to zero here if Vna is used
        if (isstress || isforce)
        {
            if (XC_Functional::get_func_type() == 3)
            {
                Gint_inout inout(DM_R,
                                 GlobalC::pot.vr_eff1,
                                 GlobalC::pot.vofk_eff1,
                                 isforce,
                                 isstress,
                                 &fvl_dphi,
                                 &svl_dphi,
                                 Gint_Tools::job_type::force_meta);
                this->UHM->GK.cal_gint(&inout);
            }
            else
            {
                Gint_inout inout(DM_R,
                                 GlobalC::pot.vr_eff1,
                                 isforce,
                                 isstress,
                                 &fvl_dphi,
                                 &svl_dphi,
                                 Gint_Tools::job_type::force);
                this->UHM->GK.cal_gint(&inout);
            }
        }
    }

    if (isstress)
    {
        StressTools::stress_fill(-1.0, GlobalC::ucell.omega, svl_dphi);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvl_dphi_k");
    return;
}

void Force_LCAO_k_tddft::calFvnlDbeta(double** dm2d,
                                const bool& isforce,
                                const bool& isstress,
                                ModuleBase::matrix& fvnl_dbeta,
                                ModuleBase::matrix& svnl_dbeta,
                                const int& vnl_method)
{
    ModuleBase::TITLE("Force_LCAO_k", "calFvnlDbeta");
    if (GlobalV::NSPIN == 4 || vnl_method == 0)
    {
        this->cal_fvnl_dbeta_k(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta);
    }
    else if (vnl_method == 1)
    {
        this->cal_fvnl_dbeta_k_new(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Force_LCAO_k", "This method has not been implemented");
    }
}


// must consider three-center H matrix.
void Force_LCAO_k_tddft::cal_fvnl_dbeta_k(double** dm2d,
                                    const bool isforce,
                                    const bool isstress,
                                    ModuleBase::matrix& fvnl_dbeta,
                                    ModuleBase::matrix& svnl_dbeta)
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_fvnl_dbeta_k");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvnl_dbeta_k");
    const Parallel_Orbitals* pv = this->ParaV;
    int iir = 0;
    ModuleBase::Vector3<double> tau1;
    ModuleBase::Vector3<double> tau2;
    ModuleBase::Vector3<double> dtau;
    ModuleBase::Vector3<double> tau0;
    ModuleBase::Vector3<double> dtau1;
    ModuleBase::Vector3<double> dtau2;

    double rcut;
    double distance;

    double rcut1;
    double rcut2;
    double distance1;
    double distance2;

    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        const Atom* atom1 = &GlobalC::ucell.atoms[T1];

        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            // GlobalC::GridD.Find_atom( tau1 );
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            // const int iat = GlobalC::ucell.itia2iat(T1, I1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

            for (int ad2 = 0; ad2 < GlobalC::GridD.getAdjacentNum() + 1; ++ad2)
            {
                const int T2 = GlobalC::GridD.getType(ad2);
                const Atom* atom2 = &GlobalC::ucell.atoms[T2];
                const int I2 = GlobalC::GridD.getNatom(ad2);
                // const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = GlobalC::GridD.getAdjacentTau(ad2);

                dtau = tau2 - tau1;
                distance = dtau.norm() * GlobalC::ucell.lat0;
                rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut();

                // check if this a adjacent atoms.
                bool is_adj = false;
                if (distance < rcut)
                    is_adj = true;
                else if (distance >= rcut)
                {
                    for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum() + 1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);
                        if (GlobalC::ucell.infoNL.nproj[T0] == 0)
                            continue;
                        const int I0 = GlobalC::GridD.getNatom(ad0);
                        // const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                        // const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        dtau2 = tau0 - tau2;
                        distance2 = dtau2.norm() * GlobalC::ucell.lat0;
                        rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if (distance1 < rcut1 && distance2 < rcut2)
                        {
                            is_adj = true;
                            break;
                        }
                    }
                }

                if (is_adj)
                {
                    // < psi1 | all projectors | psi2 >
                    // ----------------------------- enter the nnr increaing zone -------------------------
                    for (int j = 0; j < atom1->nw; ++j)
                    {
                        const int iw1_all = start1 + j;
                        const int mu = pv->trace_loc_row[iw1_all];
                        if (mu < 0)
                            continue;
                        for (int k = 0; k < atom2->nw; ++k)
                        {
                            const int iw2_all = start2 + k;
                            const int nu = pv->trace_loc_col[iw2_all];
                            if (nu < 0)
                                continue;

                            for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum() + 1; ++ad0)
                            {
                                const int T0 = GlobalC::GridD.getType(ad0);
                                if (GlobalC::ucell.infoNL.nproj[T0] == 0)
                                    continue;
                                const int I0 = GlobalC::GridD.getNatom(ad0);
                                const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                                // const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
                                tau0 = GlobalC::GridD.getAdjacentTau(ad0);

                                dtau1 = tau0 - tau1;
                                distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                                rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                                dtau2 = tau0 - tau2;
                                distance2 = dtau2.norm() * GlobalC::ucell.lat0;
                                rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                                double r0[3];
                                double r1[3];
                                r1[0] = (tau1.x - tau0.x);
                                r1[1] = (tau1.y - tau0.y);
                                r1[2] = (tau1.z - tau0.z);
                                r0[0] = (tau2.x - tau0.x);
                                r0[1] = (tau2.y - tau0.y);
                                r0[2] = (tau2.z - tau0.z);

                                if (distance1 < rcut1 && distance2 < rcut2)
                                {
                                    // const Atom* atom0 = &GlobalC::ucell.atoms[T0];
                                    double nlm[3] = {0, 0, 0};

                                    GlobalC::UOT.snap_psibeta(
                                        GlobalC::ORB,
                                        GlobalC::ucell.infoNL,
                                        nlm,
                                        1,
                                        tau2,
                                        T2,
                                        atom2->iw2l[k], // L2
                                        atom2->iw2m[k], // m2
                                        atom2->iw2n[k], // n2
                                        tau1,
                                        T1,
                                        atom1->iw2l[j], // L1
                                        atom1->iw2m[j], // m1
                                        atom1->iw2n[j], // N1
                                        tau0,
                                        T0,
                                        GlobalC::ucell.atoms[T0].dion,
                                        GlobalV::NSPIN,
                                        GlobalC::ucell.atoms[T0].d_so,
                                        GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
                                        GlobalC::ucell.atoms[T0].index1_soc[0],
                                        GlobalC::ucell.atoms[T0].index2_soc[0],
                                        GlobalC::ucell.atoms[T0].nproj_soc); // mohan  add 2021-05-07

                                    double nlm1[3] = {0, 0, 0};
                                    if (isstress)
                                    {
                                        GlobalC::UOT.snap_psibeta(
                                            GlobalC::ORB,
                                            GlobalC::ucell.infoNL,
                                            nlm1,
                                            1,
                                            tau1,
                                            T1,
                                            atom1->iw2l[j], // L1
                                            atom1->iw2m[j], // m1
                                            atom1->iw2n[j], // n1
                                            tau2,
                                            T2,
                                            atom2->iw2l[k], // L2
                                            atom2->iw2m[k], // m2
                                            atom2->iw2n[k], // N2
                                            tau0,
                                            T0,
                                            GlobalC::ucell.atoms[T0].dion,
                                            GlobalV::NSPIN,
                                            GlobalC::ucell.atoms[T0].d_so,
                                            GlobalC::ucell.atoms[T0].non_zero_count_soc[0], // index stands for spin
                                            GlobalC::ucell.atoms[T0].index1_soc[0],
                                            GlobalC::ucell.atoms[T0].index2_soc[0],
                                            GlobalC::ucell.atoms[T0].nproj_soc); // mohan  add 2021-05-07
                                    }
                                    /// only one projector for each atom force, but another projector for stress
                                    for (int is = 0; is < GlobalV::NSPIN; ++is)
                                    {
                                        double dm2d2 = 2.0 * dm2d[is][iir];
                                        for (int jpol = 0; jpol < 3; jpol++)
                                        {
                                            if (isforce)
                                            {
                                                fvnl_dbeta(iat0, jpol) -= dm2d2 * nlm[jpol];
                                            }
                                            if (isstress)
                                            {
                                                for (int ipol = jpol; ipol < 3; ipol++)
                                                {
                                                    svnl_dbeta(jpol, ipol)
                                                        += dm2d[is][iir]
                                                           * (nlm[jpol] * r1[ipol] + nlm1[jpol] * r0[ipol]);
                                                }
                                            }
                                        }
                                    }

                                } // distance
                            } // ad0

                            ++iir;
                        } // k
                    } // j
                } // distance
            } // ad2
        } // I1
    } // T1

    assert(iir == pv->nnr);

    if (isstress)
    {
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, svnl_dbeta);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvnl_dbeta_k");
    return;
}

typedef std::tuple<int, int, int, int> key_tuple;

// must consider three-center H matrix.
void Force_LCAO_k_tddft::cal_fvnl_dbeta_k_new(double** dm2d,
                                        const bool isforce,
                                        const bool isstress,
                                        ModuleBase::matrix& fvnl_dbeta,
                                        ModuleBase::matrix& svnl_dbeta)
{
    ModuleBase::TITLE("Force_LCAO_k", "cal_fvnl_dbeta_k_new");
    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvnl_dbeta_k_new");
    const Parallel_Orbitals* pv = this->ParaV;

    // Data structure for storing <psi|beta>, for a detailed description
    // check out the same data structure in build_Nonlocal_mu_new
    std::vector<std::map<key_tuple, std::unordered_map<int, std::vector<std::vector<double>>>>> nlm_tot;

    nlm_tot.resize(GlobalC::ucell.nat);

    for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
    {

        const int it = GlobalC::ucell.iat2it[iat];
        const int ia = GlobalC::ucell.iat2ia[iat];

        // Step 1 : generate <psi|beta>
        // type of atom; distance; atomic basis; projectors

        const double Rcut_Beta = GlobalC::ucell.infoNL.Beta[it].get_rcut_max();
        const ModuleBase::Vector3<double> tau = GlobalC::ucell.atoms[it].tau[ia];
        GlobalC::GridD.Find_atom(GlobalC::ucell, tau, it, ia);

        nlm_tot[iat].clear();

        for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum() + 1; ++ad)
        {
            const int T1 = GlobalC::GridD.getType(ad);
            const int I1 = GlobalC::GridD.getNatom(ad);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
            const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

            const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad);
            const Atom* atom1 = &GlobalC::ucell.atoms[T1];
            const int nw1_tot = atom1->nw * GlobalV::NPOL;

            const ModuleBase::Vector3<double> dtau = tau1 - tau;
            const double dist1 = dtau.norm2() * pow(GlobalC::ucell.lat0, 2);
            if (dist1 > pow(Rcut_Beta + Rcut_AO1, 2))
            {
                continue;
            }

            std::unordered_map<int, std::vector<std::vector<double>>> nlm_cur;
            nlm_cur.clear();

            for (int iw1 = 0; iw1 < nw1_tot; ++iw1)
            {
                const int iw1_all = start1 + iw1;
                const int iw1_local = pv->trace_loc_row[iw1_all];
                const int iw2_local = pv->trace_loc_col[iw1_all];
                if (iw1_local < 0 && iw2_local < 0)
                    continue;
                const int iw1_0 = iw1 / GlobalV::NPOL;
                std::vector<std::vector<double>> nlm;
                GlobalC::UOT.snap_psibeta_half(GlobalC::ORB,
                                               GlobalC::ucell.infoNL,
                                               nlm,
                                               tau1,
                                               T1,
                                               atom1->iw2l[iw1_0], // L1
                                               atom1->iw2m[iw1_0], // m1
                                               atom1->iw2n[iw1_0], // N1
                                               tau,
                                               it,
                                               1); // R0,T0

                nlm_cur.insert({iw1_all, nlm});
            } // end iw
            const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
            const int rx1 = GlobalC::GridD.getBox(ad).x;
            const int ry1 = GlobalC::GridD.getBox(ad).y;
            const int rz1 = GlobalC::GridD.getBox(ad).z;
            key_tuple key_1(iat1, rx1, ry1, rz1);
            nlm_tot[iat][key_1] = nlm_cur;
        } // end ad
    }

    //=======================================================
    // Step2:
    // calculate sum_(L0,M0) beta<psi_i|beta><beta|psi_j>
    // and accumulate the value to Hloc_fixedR(i,j)
    //=======================================================
    int nnr = 0;
    ModuleBase::Vector3<double> tau1;
    ModuleBase::Vector3<double> tau2;
    ModuleBase::Vector3<double> dtau;
    ModuleBase::Vector3<double> tau0;
    ModuleBase::Vector3<double> dtau1;
    ModuleBase::Vector3<double> dtau2;

    double rcut;
    double distance;

    double rcut1;
    double rcut2;
    double distance1;
    double distance2;

    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        const Atom* atom1 = &GlobalC::ucell.atoms[T1];

        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];

            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);
            const int iat1 = GlobalC::ucell.itia2iat(T1, I1);
            const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);

            for (int ad2 = 0; ad2 < GlobalC::GridD.getAdjacentNum() + 1; ++ad2)
            {
                const int T2 = GlobalC::GridD.getType(ad2);
                const Atom* atom2 = &GlobalC::ucell.atoms[T2];

                const int I2 = GlobalC::GridD.getNatom(ad2);
                const int iat2 = GlobalC::ucell.itia2iat(T2, I2);
                const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
                tau2 = GlobalC::GridD.getAdjacentTau(ad2);

                const int rx2 = GlobalC::GridD.getBox(ad2).x;
                const int ry2 = GlobalC::GridD.getBox(ad2).y;
                const int rz2 = GlobalC::GridD.getBox(ad2).z;

                dtau = tau2 - tau1;
                distance = dtau.norm2() * pow(GlobalC::ucell.lat0, 2);
                // this rcut is in order to make nnr consistent
                // with other matrix.
                rcut = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Phi[T2].getRcut(), 2);

                // check if this a adjacent atoms.
                bool is_adj = false;
                if (distance < rcut)
                    is_adj = true;
                else if (distance >= rcut)
                {
                    for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum() + 1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);
                        if (GlobalC::ucell.infoNL.nproj[T0] == 0)
                            continue;
                        const int I0 = GlobalC::GridD.getNatom(ad0);
                        // const int iat0 = GlobalC::ucell.itia2iat(T0, I0);
                        // const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);

                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);
                        dtau1 = tau0 - tau1;
                        distance1 = dtau1.norm() * GlobalC::ucell.lat0;
                        rcut1 = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        dtau2 = tau0 - tau2;
                        distance2 = dtau2.norm() * GlobalC::ucell.lat0;
                        rcut2 = GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max();

                        if (distance1 < rcut1 && distance2 < rcut2)
                        {
                            is_adj = true;
                            break;
                        }
                    }
                }

                if (is_adj)
                {
                    for (int ad0 = 0; ad0 < GlobalC::GridD.getAdjacentNum() + 1; ++ad0)
                    {
                        const int T0 = GlobalC::GridD.getType(ad0);
                        const int I0 = GlobalC::GridD.getNatom(ad0);
                        const int iat = GlobalC::ucell.itia2iat(T0, I0);

                        // mohan add 2010-12-19
                        if (GlobalC::ucell.infoNL.nproj[T0] == 0)
                            continue;

                        // const int I0 = GlobalC::GridD.getNatom(ad0);
                        // const int start0 = GlobalC::ucell.itiaiw2iwt(T0, I0, 0);
                        tau0 = GlobalC::GridD.getAdjacentTau(ad0);

                        dtau1 = tau0 - tau1;
                        dtau2 = tau0 - tau2;
                        const double distance1 = dtau1.norm2() * pow(GlobalC::ucell.lat0, 2);
                        const double distance2 = dtau2.norm2() * pow(GlobalC::ucell.lat0, 2);

                        // seems a bug here!! mohan 2011-06-17
                        rcut1 = pow(GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(), 2);
                        rcut2 = pow(GlobalC::ORB.Phi[T2].getRcut() + GlobalC::ucell.infoNL.Beta[T0].get_rcut_max(), 2);

                        double r0[3];
                        double r1[3];
                        r1[0] = (tau1.x - tau0.x);
                        r1[1] = (tau1.y - tau0.y);
                        r1[2] = (tau1.z - tau0.z);
                        r0[0] = (tau2.x - tau0.x);
                        r0[1] = (tau2.y - tau0.y);
                        r0[2] = (tau2.z - tau0.z);

                        if (distance1 >= rcut1 || distance2 >= rcut2)
                        {
                            continue;
                        }

                        const int rx0 = GlobalC::GridD.getBox(ad0).x;
                        const int ry0 = GlobalC::GridD.getBox(ad0).y;
                        const int rz0 = GlobalC::GridD.getBox(ad0).z;
                        key_tuple key1(iat1, -rx0, -ry0, -rz0);
                        key_tuple key2(iat2, rx2 - rx0, ry2 - ry0, rz2 - rz0);

                        int nnr_inner = 0;
                        for (int j = 0; j < atom1->nw * GlobalV::NPOL; j++)
                        {
                            const int j0 = j / GlobalV::NPOL; // added by zhengdy-soc
                            const int iw1_all = start1 + j;
                            const int mu = pv->trace_loc_row[iw1_all];
                            if (mu < 0)
                                continue;

                            for (int k = 0; k < atom2->nw * GlobalV::NPOL; k++)
                            {
                                const int k0 = k / GlobalV::NPOL;
                                const int iw2_all = start2 + k;
                                const int nu = pv->trace_loc_col[iw2_all];
                                if (nu < 0)
                                    continue;

                                // const Atom* atom0 = &GlobalC::ucell.atoms[T0];
                                double nlm[3] = {0, 0, 0};
                                std::vector<double> nlm_1 = nlm_tot[iat][key2][iw2_all][0];
                                std::vector<std::vector<double>> nlm_2;
                                nlm_2.resize(3);
                                for (int i = 0; i < 3; i++)
                                {
                                    nlm_2[i] = nlm_tot[iat][key1][iw1_all][i + 1];
                                }

                                assert(nlm_1.size() == nlm_2[0].size());

                                const int nproj = GlobalC::ucell.infoNL.nproj[T0];
                                int ib = 0;
                                for (int nb = 0; nb < nproj; nb++)
                                {
                                    const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
                                    for (int m = 0; m < 2 * L0 + 1; m++)
                                    {
                                        for (int ir = 0; ir < 3; ir++)
                                        {
                                            nlm[ir]
                                                += nlm_2[ir][ib] * nlm_1[ib] * GlobalC::ucell.atoms[T0].dion(nb, nb);
                                        }
                                        ib += 1;
                                    }
                                }
                                assert(ib == nlm_1.size());

                                double nlm1[3] = {0, 0, 0};
                                if (isstress)
                                {
                                    std::vector<double> nlm_1 = nlm_tot[iat][key1][iw1_all][0];
                                    std::vector<std::vector<double>> nlm_2;
                                    nlm_2.resize(3);
                                    for (int i = 0; i < 3; i++)
                                    {
                                        nlm_2[i] = nlm_tot[iat][key2][iw2_all][i + 1];
                                    }

                                    assert(nlm_1.size() == nlm_2[0].size());

                                    const int nproj = GlobalC::ucell.infoNL.nproj[T0];
                                    int ib = 0;
                                    for (int nb = 0; nb < nproj; nb++)
                                    {
                                        const int L0 = GlobalC::ucell.infoNL.Beta[T0].Proj[nb].getL();
                                        for (int m = 0; m < 2 * L0 + 1; m++)
                                        {
                                            for (int ir = 0; ir < 3; ir++)
                                            {
                                                nlm1[ir] += nlm_2[ir][ib] * nlm_1[ib]
                                                            * GlobalC::ucell.atoms[T0].dion(nb, nb);
                                            }
                                            ib += 1;
                                        }
                                    }
                                    assert(ib == nlm_1.size());
                                }

                                /// only one projector for each atom force, but another projector for stress
                                for (int is = 0; is < GlobalV::NSPIN; ++is)
                                {
                                    double dm2d2 = 2.0 * dm2d[is][nnr + nnr_inner];
                                    for (int jpol = 0; jpol < 3; jpol++)
                                    {
                                        if (isforce)
                                        {
                                            fvnl_dbeta(iat, jpol) -= dm2d2 * nlm[jpol];
                                        }
                                        if (isstress)
                                        {
                                            for (int ipol = jpol; ipol < 3; ipol++)
                                            {
                                                svnl_dbeta(jpol, ipol)
                                                    += dm2d[is][nnr + nnr_inner]
                                                       * (nlm[jpol] * r1[ipol] + nlm1[jpol] * r0[ipol]);
                                            }
                                        }
                                    }
                                }
                                nnr_inner++;
                            } // k
                        } // j
                    } // ad0

                    // outer circle : accumulate nnr
                    for (int j = 0; j < atom1->nw * GlobalV::NPOL; j++)
                    {
                        const int j0 = j / GlobalV::NPOL; // added by zhengdy-soc
                        const int iw1_all = start1 + j;
                        const int mu = pv->trace_loc_row[iw1_all];
                        if (mu < 0)
                            continue;

                        // fix a serious bug: atom2[T2] -> atom2
                        // mohan 2010-12-20
                        for (int k = 0; k < atom2->nw * GlobalV::NPOL; k++)
                        {
                            const int k0 = k / GlobalV::NPOL;
                            const int iw2_all = start2 + k;
                            const int nu = pv->trace_loc_col[iw2_all];
                            if (nu < 0)
                                continue;

                            nnr++;
                        }
                    }
                } // is_adj
            } // ad2
        } // I1
    } // T1

    assert(nnr == pv->nnr);

    if (isstress)
    {
        StressTools::stress_fill(GlobalC::ucell.lat0, GlobalC::ucell.omega, svnl_dbeta);
    }

    ModuleBase::timer::tick("Force_LCAO_k", "cal_fvnl_dbeta_k_new");
    return;
}*/