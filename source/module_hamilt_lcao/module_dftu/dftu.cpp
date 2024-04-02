//==========================================================
// Author:Xin Qu
// DATE : 2019-12-10
//==========================================================
#include "dftu.h"

#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/inverse_matrix.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/magnetism.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_base/scalapack_connector.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>

namespace GlobalC
{
ModuleDFTU::DFTU dftu;
}

namespace ModuleDFTU
{
DFTU::DFTU()
{
}

DFTU::~DFTU()
{
}

void DFTU::init(UnitCell& cell, // unitcell class
                LCAO_Matrix& lm,
                const int& nks)
{
    ModuleBase::TITLE("DFTU", "init");

#ifndef __MPI
    std::cout << "DFT+U module is only accessible in mpi versioin" << std::endl;
    exit(0);
#endif

    this->LM = &lm;

    // needs reconstructions in future
    // global parameters, need to be removed in future
    const int npol = GlobalV::NPOL; // number of polarization directions
    const int nlocal = GlobalV::NLOCAL; // number of total local orbitals
    const int nspin = GlobalV::NSPIN; // number of spins

    /*
    not implemented yet, no need to check for now
    if(dftu_type==1 && double_counting==1) cal_type = 1;
    else if(dftu_type==1 && double_counting==2) cal_type = 2;
    else if(dftu_type==2 && double_counting==1) cal_type = 3;
    else if(dftu_type==2 && double_counting==2) cal_type = 4;
    else ModuleBase::WARNING_QUIT("DFT+U", "Wrong parameter");

    if(cal_type!=3)  ModuleBase::WARNING_QUIT("DFT+U", "Not available yet!");
    */

    this->EU = 0.0;

    this->locale.resize(cell.nat);
    this->locale_save.resize(cell.nat);

    this->iatlnmipol2iwt.resize(cell.nat);

    int num_locale = 0;
    // it:index of type of atom
    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            // ia:index of atoms of this type
            // determine the size of locale
            const int iat = cell.itia2iat(it, ia);

            locale[iat].resize(cell.atoms[it].nwl + 1);
            locale_save[iat].resize(cell.atoms[it].nwl + 1);

            for (int l = 0; l <= cell.atoms[it].nwl; l++)
            {
                const int N = cell.atoms[it].l_nchi[l];

                locale[iat][l].resize(N);
                locale_save[iat][l].resize(N);

                for (int n = 0; n < N; n++)
                {
                    if (nspin == 1 || nspin == 2)
                    {
                        locale[iat][l][n].resize(2);
                        locale_save[iat][l][n].resize(2);

                        locale[iat][l][n][0].create(2 * l + 1, 2 * l + 1);
                        locale[iat][l][n][1].create(2 * l + 1, 2 * l + 1);

                        locale_save[iat][l][n][0].create(2 * l + 1, 2 * l + 1);
                        locale_save[iat][l][n][1].create(2 * l + 1, 2 * l + 1);
                        num_locale += (2 * l + 1) * (2 * l + 1) * 2;
                    }
                    else if (nspin == 4) // SOC
                    {
                        locale[iat][l][n].resize(1);
                        locale_save[iat][l][n].resize(1);

                        locale[iat][l][n][0].create((2 * l + 1) * npol, (2 * l + 1) * npol);
                        locale_save[iat][l][n][0].create((2 * l + 1) * npol, (2 * l + 1) * npol);
                        num_locale += (2 * l + 1) * (2 * l + 1) * npol * npol;
                    }
                }
            }

            // initialize the arrry iatlnm2iwt[iat][l][n][m]
            this->iatlnmipol2iwt[iat].resize(cell.atoms[it].nwl + 1);
            for (int L = 0; L <= cell.atoms[it].nwl; L++)
            {
                this->iatlnmipol2iwt[iat][L].resize(cell.atoms[it].l_nchi[L]);

                for (int n = 0; n < cell.atoms[it].l_nchi[L]; n++)
                {
                    this->iatlnmipol2iwt[iat][L][n].resize(2 * L + 1);

                    for (int m = 0; m < 2 * L + 1; m++)
                    {
                        this->iatlnmipol2iwt[iat][L][n][m].resize(npol);
                    }
                }
            }

            for (int iw = 0; iw < cell.atoms[it].nw * npol; iw++)
            {
                int iw0 = iw / npol;
                int ipol = iw % npol;
                int iwt = cell.itiaiw2iwt(it, ia, iw);
                int l = cell.atoms[it].iw2l[iw0];
                int n = cell.atoms[it].iw2n[iw0];
                int m = cell.atoms[it].iw2m[iw0];

                this->iatlnmipol2iwt[iat][l][n][m][ipol] = iwt;
            }
        }
    }

    if (Yukawa)
    {
        this->Fk.resize(cell.ntype);

        this->U_Yukawa.resize(cell.ntype);
        this->J_Yukawa.resize(cell.ntype);

        for (int it = 0; it < cell.ntype; it++)
        {
            const int NL = cell.atoms[it].nwl + 1;

            this->Fk[it].resize(NL);
            this->U_Yukawa[it].resize(NL);
            this->J_Yukawa[it].resize(NL);

            for (int l = 0; l < NL; l++)
            {
                int N = cell.atoms[it].l_nchi[l];

                this->Fk[it][l].resize(N);
                for (int n = 0; n < N; n++)
                {
                    this->Fk[it][l][n].resize(l + 1, 0.0);
                }

                this->U_Yukawa[it][l].resize(N, 0.0);
                this->J_Yukawa[it][l].resize(N, 0.0);
            }
        }
    }

    if (omc != 0)
    {
        std::stringstream sst;
        sst << "initial_onsite.dm";
        this->read_occup_m(sst.str());
#ifdef __MPI
        this->local_occup_bcast();
#endif

        initialed_locale = true;
        this->copy_locale();
    }
    else
    {
        if (GlobalV::init_chg == "file")
        {
            std::stringstream sst;
            sst << GlobalV::global_out_dir << "onsite.dm";
            this->read_occup_m(sst.str());
#ifdef __MPI
            this->local_occup_bcast();
#endif
            initialed_locale = true;
        }
        else
        {
            this->zero_locale();
        }
    }

    ModuleBase::Memory::record("DFTU::locale", sizeof(double) * num_locale);
    return;
}

void DFTU::cal_energy_correction(const int istep)
{
    ModuleBase::TITLE("DFTU", "cal_energy_correction");
    ModuleBase::timer::tick("DFTU", "cal_energy_correction");
    if (!initialed_locale)
    {
        ModuleBase::timer::tick("DFTU", "cal_energy_correction");
        return;
    }
    this->EU = 0.0;
    double EU_dc = 0.0;

    for (int T = 0; T < GlobalC::ucell.ntype; T++)
    {
        const int NL = GlobalC::ucell.atoms[T].nwl + 1;
        const int LC = orbital_corr[T];
        for (int I = 0; I < GlobalC::ucell.atoms[T].na; I++)
        {
            if (LC == -1)
                continue;

            const int iat = GlobalC::ucell.itia2iat(T, I);
            const int L = orbital_corr[T];

            for (int l = 0; l < NL; l++)
            {
                if (l != orbital_corr[T])
                    continue;

                const int N = GlobalC::ucell.atoms[T].l_nchi[l];

                const int m_tot = 2 * l + 1;

                // part 1: calculate the DFT+U energy correction
                for (int n = 0; n < N; n++)
                {
                    if (n != 0)
                        continue;

                    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
                    {
                        for (int spin = 0; spin < 2; spin++)
                        {
                            double nm_trace = 0.0;
                            double nm2_trace = 0.0;

                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                nm_trace += this->locale[iat][l][n][spin](m0, m0);
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    nm2_trace += this->locale[iat][l][n][spin](m0, m1)
                                                 * this->locale[iat][l][n][spin](m1, m0);
                                }
                            }
                            if (Yukawa)
                                this->EU += 0.5 * (this->U_Yukawa[T][l][n] - this->J_Yukawa[T][l][n])
                                            * (nm_trace - nm2_trace);
                            else
                                this->EU += 0.5 * this->U[T] * (nm_trace - nm2_trace);
                        }
                    }
                    else if (GlobalV::NSPIN == 4) // SOC
                    {
                        double nm_trace = 0.0;
                        double nm2_trace = 0.0;

                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < GlobalV::NPOL; ipol0++)
                            {
                                const int m0_all = m0 + (2 * l + 1) * ipol0;
                                nm_trace += this->locale[iat][l][n][0](m0_all, m0_all);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < GlobalV::NPOL; ipol1++)
                                    {
                                        int m1_all = m1 + (2 * l + 1) * ipol1;

                                        nm2_trace += this->locale[iat][l][n][0](m0_all, m1_all)
                                                     * this->locale[iat][l][n][0](m1_all, m0_all);
                                    }
                                }
                            }
                        }
                        if (Yukawa)
                            this->EU
                                += 0.5 * (this->U_Yukawa[T][l][n] - this->J_Yukawa[T][l][n]) * (nm_trace - nm2_trace);
                        else
                            this->EU += 0.5 * this->U[T] * (nm_trace - nm2_trace);
                    }

                    // calculate the double counting term included in eband
                    for (int m1 = 0; m1 < 2 * l + 1; m1++)
                    {
                        for (int ipol1 = 0; ipol1 < GlobalV::NPOL; ipol1++)
                        {
                            const int m1_all = m1 + ipol1 * (2 * l + 1);
                            for (int m2 = 0; m2 < 2 * l + 1; m2++)
                            {
                                for (int ipol2 = 0; ipol2 < GlobalV::NPOL; ipol2++)
                                {
                                    const int m2_all = m2 + ipol2 * (2 * l + 1);

                                    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
                                    {
                                        for (int is = 0; is < 2; is++)
                                        {
                                            double VU = 0.0;
                                            VU = get_onebody_eff_pot(T, iat, l, n, is, m1_all, m2_all, 0);
                                            EU_dc += VU * this->locale[iat][l][n][is](m1_all, m2_all);
                                        }
                                    }
                                    else if (GlobalV::NSPIN == 4) // SOC
                                    {
                                        double VU = 0.0;
                                        VU = get_onebody_eff_pot(T, iat, l, n, 0, m1_all, m2_all, 0);
                                        EU_dc += VU * this->locale[iat][l][n][0](m1_all, m2_all);
                                    }
                                }
                            }
                        }
                    }
                } // end n
            } // end L
        } // end I
    } // end T

    // substract the double counting EU_dc included in band energy eband
    this->EU -= EU_dc;

    ModuleBase::timer::tick("DFTU", "cal_energy_correction");
    return;
}

void DFTU::uramping_update()
{
    // if uramping < 0.1, use the original U
    if(this->uramping < 0.1) return;
    // loop to change U
    for(int i = 0; i < this->U0.size(); i++)
    {
        if (this->U[i] + this->uramping < this->U0[i] ) 
        {
            this->U[i] += this->uramping;
        }
        else
        {
            this->U[i] = this->U0[i];
        }
    }
}

bool DFTU::u_converged()
{
    for(int i = 0; i < this->U0.size(); i++)
    {
        if (this->U[i] != this->U0[i]) 
        {
            return false;
        }
    }
    return true;
}

} // namespace ModuleDFTU
