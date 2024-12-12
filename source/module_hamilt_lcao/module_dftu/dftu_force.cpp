//==========================================================
// Author:Xin Qu
#include "module_parameter/parameter.h"
// DATE : 2019-12-10
//==========================================================
#ifdef __LCAO
#include "dftu.h"
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/inverse_matrix.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/magnetism.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>

extern "C"
{
    // I'm not sure what's happenig here, but the interface in scalapack_connecter.h
    // does not seem to work, so I'll use this one here
    void pzgemm_(const char* transa,
                 const char* transb,
                 const int* M,
                 const int* N,
                 const int* K,
                 const std::complex<double>* alpha,
                 const std::complex<double>* A,
                 const int* IA,
                 const int* JA,
                 const int* DESCA,
                 const std::complex<double>* B,
                 const int* IB,
                 const int* JB,
                 const int* DESCB,
                 const std::complex<double>* beta,
                 std::complex<double>* C,
                 const int* IC,
                 const int* JC,
                 const int* DESCC);

    void pdgemm_(const char* transa,
                 const char* transb,
                 const int* M,
                 const int* N,
                 const int* K,
                 const double* alpha,
                 const double* A,
                 const int* IA,
                 const int* JA,
                 const int* DESCA,
                 const double* B,
                 const int* IB,
                 const int* JB,
                 const int* DESCB,
                 const double* beta,
                 double* C,
                 const int* IC,
                 const int* JC,
                 const int* DESCC);
}

namespace ModuleDFTU
{

void DFTU::force_stress(const UnitCell& ucell,
                        const Grid_Driver& gd,
                        const elecstate::ElecState* pelec,
                        const Parallel_Orbitals& pv,
                        ForceStressArrays& fsr, // mohan add 2024-06-16
                        ModuleBase::matrix& force_dftu,
                        ModuleBase::matrix& stress_dftu,
                        const K_Vectors& kv)
{
    ModuleBase::TITLE("DFTU", "force_stress");
    ModuleBase::timer::tick("DFTU", "force_stress");

    const int nlocal = PARAM.globalv.nlocal;

    if (PARAM.inp.cal_force)
    {
        force_dftu.zero_out();
    }
    if (PARAM.inp.cal_stress)
    {
        stress_dftu.zero_out();
    }

    if (PARAM.globalv.gamma_only_local)
    {
        const char transN = 'N';
        const char transT = 'T';
        const int one_int = 1;
        const double alpha = 1.0;
        const double beta = 0.0;

        std::vector<double> rho_VU(pv.nloc);

        for (int ik = 0; ik < kv.get_nks(); ik++)
        {

            const int spin = kv.isk[ik];

            double* VU = new double[pv.nloc];

            this->cal_VU_pot_mat_real(spin, false, VU);

            const std::vector<std::vector<double>>& dmk
                = dynamic_cast<const elecstate::ElecStateLCAO<double>*>(pelec)->get_DM()->get_DMK_vector();

#ifdef __MPI
            pdgemm_(&transT,
                    &transN,
                    &nlocal,
                    &nlocal,
                    &nlocal,
                    &alpha,
                    dmk[spin].data(),
                    &one_int,
                    &one_int,
                    pv.desc,
                    VU,
                    &one_int,
                    &one_int,
                    pv.desc,
                    &beta,
                    &rho_VU[0],
                    &one_int,
                    &one_int,
                    pv.desc);
#endif

            delete[] VU;

            if (PARAM.inp.cal_force)
            {
                this->cal_force_gamma(ucell,&rho_VU[0], pv, fsr.DSloc_x, fsr.DSloc_y, fsr.DSloc_z, force_dftu);
            }

            if (PARAM.inp.cal_stress)
            {
                this->cal_stress_gamma(ucell,
                                       pv,
                                       &gd,
                                       fsr.DSloc_x,
                                       fsr.DSloc_y,
                                       fsr.DSloc_z,
                                       fsr.DH_r,
                                       &rho_VU[0],
                                       stress_dftu);
            }
        } // ik
    }
    else
    {
        const char transN = 'N';
        const char transT = 'T';
        const int one_int = 1;
        const std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);

        std::vector<std::complex<double>> rho_VU(pv.nloc);

        for (int ik = 0; ik < kv.get_nks(); ik++)
        {
            const int spin = kv.isk[ik];

            std::complex<double>* VU = new std::complex<double>[pv.nloc];

            this->cal_VU_pot_mat_complex(spin, false, VU);

            const std::vector<std::vector<std::complex<double>>>& dmk
                = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(pelec)
                      ->get_DM()
                      ->get_DMK_vector();

#ifdef __MPI
            pzgemm_(&transT,
                    &transN,
                    &nlocal,
                    &nlocal,
                    &nlocal,
                    &alpha,
                    dmk[ik].data(),
                    &one_int,
                    &one_int,
                    pv.desc,
                    VU,
                    &one_int,
                    &one_int,
                    pv.desc,
                    &beta,
                    &rho_VU[0],
                    &one_int,
                    &one_int,
                    pv.desc);
#endif

            delete[] VU;

            if (PARAM.inp.cal_force)
            {
                cal_force_k(ucell, gd, fsr, pv, ik, &rho_VU[0], force_dftu, kv.kvec_d);
            }
            if (PARAM.inp.cal_stress)
            {
                cal_stress_k(ucell, gd, fsr, pv, ik, &rho_VU[0], stress_dftu, kv.kvec_d);
            }
        } // ik
    }

    if (PARAM.inp.cal_force)
    {
        Parallel_Reduce::reduce_pool(force_dftu.c, force_dftu.nr * force_dftu.nc);
    }

    if (PARAM.inp.cal_stress)
    {
        Parallel_Reduce::reduce_pool(stress_dftu.c, stress_dftu.nr * stress_dftu.nc);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (i > j)
                    stress_dftu(i, j) = stress_dftu(j, i);
            }
        }

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                stress_dftu(i, j) *= ucell.lat0 / ucell.omega;
            }
        }
    }
    ModuleBase::timer::tick("DFTU", "force_stress");

    return;
}

void DFTU::cal_force_k(const UnitCell& ucell,
                       const Grid_Driver& gd,
                       ForceStressArrays& fsr,
                       const Parallel_Orbitals& pv,
                       const int ik,
                       const std::complex<double>* rho_VU,
                       ModuleBase::matrix& force_dftu,
                       const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("DFTU", "cal_force_k");
    ModuleBase::timer::tick("DFTU", "cal_force_k");

    const char transN = 'N';
    const char transC = 'C';
    const int one_int = 1;
    const std::complex<double> zero(0.0, 0.0);
    const std::complex<double> one(1.0, 0.0);

    std::vector<std::complex<double>> dm_VU_dSm(pv.nloc);
    std::vector<std::complex<double>> dSm_k(pv.nloc);

    for (int dim = 0; dim < 3; dim++)
    {
        this->folding_matrix_k(ucell, gd, fsr, pv, ik, dim + 1, 0, &dSm_k[0], kvec_d);

#ifdef __MPI
        pzgemm_(&transN,
                &transC,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &one,
                &dSm_k[0],
                &one_int,
                &one_int,
                pv.desc,
                rho_VU,
                &one_int,
                &one_int,
                pv.desc,
                &zero,
                &dm_VU_dSm[0],
                &one_int,
                &one_int,
                pv.desc);
#endif

        for (int ir = 0; ir < pv.nrow; ir++)
        {
            const int iwt1 = pv.local2global_row(ir);
            const int iat1 = ucell.iwt2iat[iwt1];

            for (int ic = 0; ic < pv.ncol; ic++)
            {
                const int iwt2 = pv.local2global_col(ic);
                const int irc = ic * pv.nrow + ir;

                if (iwt1 == iwt2)
                    force_dftu(iat1, dim) += dm_VU_dSm[irc].real();

            } // end ic
        }     // end ir

#ifdef __MPI
        pzgemm_(&transN,
                &transN,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &one,
                &dSm_k[0],
                &one_int,
                &one_int,
                pv.desc,
                rho_VU,
                &one_int,
                &one_int,
                pv.desc,
                &zero,
                &dm_VU_dSm[0],
                &one_int,
                &one_int,
                pv.desc);
#endif

        for (int it = 0; it < ucell.ntype; it++)
        {
            const int NL = ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

            if (LC == -1)
                continue;
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it])
                        continue;
                    const int N = ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0)
                            continue;

                        for (int m = 0; m < 2 * l + 1; m++)
                        {
                            for (int ipol = 0; ipol < PARAM.globalv.npol; ipol++)
                            {
                                const int iwt = this->iatlnmipol2iwt[iat][l][n][m][ipol];
                                const int mu = pv.global2local_row(iwt);
                                const int nu = pv.global2local_col(iwt);
                                if (mu < 0 || nu < 0)
                                    continue;

                                force_dftu(iat, dim) += dm_VU_dSm[nu * pv.nrow + mu].real();
                            }
                        } //
                    }     // n
                }         // l
            }             // ia
        }                 // it
    }                     // end dim
    ModuleBase::timer::tick("DFTU", "cal_force_k");

    return;
}

void DFTU::cal_stress_k(const UnitCell& ucell,
                        const Grid_Driver& gd,
                        ForceStressArrays& fsr,
                        const Parallel_Orbitals& pv,
                        const int ik,
                        const std::complex<double>* rho_VU,
                        ModuleBase::matrix& stress_dftu,
                        const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("DFTU", "cal_stress_k");
    ModuleBase::timer::tick("DFTU", "cal_stress_k");

    const int nlocal = PARAM.globalv.nlocal;

    const char transN = 'N';
    const int one_int = 1;
    const std::complex<double> minus_half(-0.5, 0.0);
    const std::complex<double> zero(0.0, 0.0);
    const std::complex<double> one(1.0, 0.0);

    std::vector<std::complex<double>> dm_VU_sover(pv.nloc);
    std::vector<std::complex<double>> dSR_k(pv.nloc);

    for (int dim1 = 0; dim1 < 3; dim1++)
    {
        for (int dim2 = dim1; dim2 < 3; dim2++)
        {
            this->folding_matrix_k(ucell, gd, fsr, pv, ik, dim1 + 4, dim2, &dSR_k[0], kvec_d);

#ifdef __MPI
            pzgemm_(&transN,
                    &transN,
                    &nlocal,
                    &nlocal,
                    &nlocal,
                    &minus_half,
                    rho_VU,
                    &one_int,
                    &one_int,
                    pv.desc,
                    &dSR_k[0],
                    &one_int,
                    &one_int,
                    pv.desc,
                    &zero,
                    &dm_VU_sover[0],
                    &one_int,
                    &one_int,
                    pv.desc);
#endif

            for (int ir = 0; ir < pv.nrow; ir++)
            {
                const int iwt1 = pv.local2global_row(ir);
                for (int ic = 0; ic < pv.ncol; ic++)
                {
                    const int iwt2 = pv.local2global_col(ic);
                    const int irc = ic * pv.nrow + ir;

                    if (iwt1 == iwt2)
                        stress_dftu(dim1, dim2) += 2.0 * dm_VU_sover[irc].real();
                } // end ic
            }     // end ir

        } // end dim2
    }     // end dim1
    ModuleBase::timer::tick("DFTU", "cal_stress_k");

    return;
}

void DFTU::cal_force_gamma(const UnitCell& ucell,
                           const double* rho_VU,
                           const Parallel_Orbitals& pv,
                           double* dsloc_x,
                           double* dsloc_y,
                           double* dsloc_z,
                           ModuleBase::matrix& force_dftu)
{
    ModuleBase::TITLE("DFTU", "cal_force_gamma");
    ModuleBase::timer::tick("DFTU", "cal_force_gamma");
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const double one = 1.0, zero = 0.0, minus_one = -1.0;

    std::vector<double> dm_VU_dSm(pv.nloc);

    for (int dim = 0; dim < 3; dim++)
    {
        double* tmp_ptr;
        if (dim == 0)
        {
            tmp_ptr = dsloc_x;
        }
        else if (dim == 1)
        {
            tmp_ptr = dsloc_y;
        }
        else if (dim == 2)
        {
            tmp_ptr = dsloc_z;
        }

#ifdef __MPI
        pdgemm_(&transN,
                &transT,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &one,
                tmp_ptr,
                &one_int,
                &one_int,
                pv.desc,
                rho_VU,
                &one_int,
                &one_int,
                pv.desc,
                &zero,
                &dm_VU_dSm[0],
                &one_int,
                &one_int,
                pv.desc);
#endif

        for (int ir = 0; ir < pv.nrow; ir++)
        {
            const int iwt1 = pv.local2global_row(ir);
            const int iat1 = ucell.iwt2iat[iwt1];

            for (int ic = 0; ic < pv.ncol; ic++)
            {
                const int iwt2 = pv.local2global_col(ic);
                const int irc = ic * pv.nrow + ir;

                if (iwt1 == iwt2)
                    force_dftu(iat1, dim) += dm_VU_dSm[irc];

            } // end ic
        }     // end ir

#ifdef __MPI
        pdgemm_(&transN,
                &transT,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &PARAM.globalv.nlocal,
                &one,
                tmp_ptr,
                &one_int,
                &one_int,
                pv.desc,
                rho_VU,
                &one_int,
                &one_int,
                pv.desc,
                &zero,
                &dm_VU_dSm[0],
                &one_int,
                &one_int,
                pv.desc);
#endif

        for (int it = 0; it < ucell.ntype; it++)
        {
            const int NL = ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

            if (LC == -1)
                continue;
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it])
                        continue;

                    const int N = ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0)
                            continue;

                        // Calculate the local occupation number matrix
                        for (int m = 0; m < 2 * l + 1; m++)
                        {
                            for (int ipol = 0; ipol < PARAM.globalv.npol; ipol++)
                            {
                                const int iwt = this->iatlnmipol2iwt[iat][l][n][m][ipol];
                                const int mu = pv.global2local_row(iwt);
                                const int nu = pv.global2local_col(iwt);
                                if (mu < 0 || nu < 0)
                                    continue;

                                force_dftu(iat, dim) += dm_VU_dSm[nu * pv.nrow + mu];
                            }
                        } //
                    }     // n
                }         // l
            }             // ia
        }                 // it

    } // end dim
    ModuleBase::timer::tick("DFTU", "cal_force_gamma");

    return;
}

void DFTU::cal_stress_gamma(const UnitCell& ucell,
                            const Parallel_Orbitals& pv,
                            const Grid_Driver* gd,
                            double* dsloc_x,
                            double* dsloc_y,
                            double* dsloc_z,
                            double* dh_r,
                            const double* rho_VU,
                            ModuleBase::matrix& stress_dftu)
{
    ModuleBase::TITLE("DFTU", "cal_stress_gamma");
    ModuleBase::timer::tick("DFTU", "cal_stress_gamma");

    const char transN = 'N';
    const int one_int = 1;
    const double zero = 0.0;
    const double minus_half = -0.5;
    const double one = 1.0;

    std::vector<double> dSR_gamma(pv.nloc);
    std::vector<double> dm_VU_sover(pv.nloc);

    const int nlocal = PARAM.globalv.nlocal;

    for (int dim1 = 0; dim1 < 3; dim1++)
    {
        for (int dim2 = dim1; dim2 < 3; dim2++)
        {
            this->fold_dSR_gamma(ucell, pv, gd, dsloc_x, dsloc_y, dsloc_z, dh_r, dim1, dim2, &dSR_gamma[0]);

#ifdef __MPI
            pdgemm_(&transN,
                    &transN,
                    &nlocal,
                    &nlocal,
                    &nlocal,
                    &minus_half,
                    rho_VU,
                    &one_int,
                    &one_int,
                    pv.desc,
                    &dSR_gamma[0],
                    &one_int,
                    &one_int,
                    pv.desc,
                    &zero,
                    &dm_VU_sover[0],
                    &one_int,
                    &one_int,
                    pv.desc);
#endif

            for (int ir = 0; ir < this->paraV->nrow; ir++)
            {
                const int iwt1 = this->paraV->local2global_row(ir);

                for (int ic = 0; ic < this->paraV->ncol; ic++)
                {
                    const int iwt2 = this->paraV->local2global_col(ic);
                    const int irc = ic * this->paraV->nrow + ir;

                    if (iwt1 == iwt2)
                        stress_dftu(dim1, dim2) += 2.0 * dm_VU_sover[irc];
                } // end ic
            }     // end ir

        } // end dim2
    }     // end dim1
    ModuleBase::timer::tick("DFTU", "cal_stress_gamma");
    return;
}
} // namespace ModuleDFTU
#endif