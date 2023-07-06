//==========================================================
// Author:Xin Qu
// DATE : 2019-12-10
//==========================================================
#include <stdio.h>
#include <string.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "dftu.h"
#include "module_base/constants.h"
#include "module_base/global_function.h"
#include "module_base/inverse_matrix.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/ORB_gen_tables.h"
#include "module_elecstate/magnetism.h"
#include "module_elecstate/module_charge/charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_matrix.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

extern "C"
{
  //I'm not sure what's happenig here, but the interface in scalapack_connecter.h
  //does not seem to work, so I'll use this one here
  void pzgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const std::complex<double> *alpha,
		const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
		const std::complex<double> *B, const int *IB, const int *JB, const int *DESCB,
		const std::complex<double> *beta,
		std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);
  
  void pdgemm_(
		const char *transa, const char *transb,
		const int *M, const int *N, const int *K,
		const double *alpha,
		const double *A, const int *IA, const int *JA, const int *DESCA,
		const double *B, const int *IB, const int *JB, const int *DESCB,
		const double *beta,
		double *C, const int *IC, const int *JC, const int *DESCC);
}

namespace ModuleDFTU
{

void DFTU::force_stress(std::vector<ModuleBase::matrix>& dm_gamma,
                        std::vector<ModuleBase::ComplexMatrix>& dm_k,
                        LCAO_Matrix& lm,
                        ModuleBase::matrix& force_dftu,
                        ModuleBase::matrix& stress_dftu,
                        const K_Vectors& kv)
{
    ModuleBase::TITLE("DFTU", "force_stress");
    ModuleBase::timer::tick("DFTU", "force_stress");

    this->LM = &lm;

    if (GlobalV::CAL_FORCE)
    {
        force_dftu.zero_out();
    }
    if (GlobalV::CAL_STRESS)
    {
        stress_dftu.zero_out();
    }

    if (GlobalV::GAMMA_ONLY_LOCAL)
    {
        const char transN = 'N', transT = 'T';
        const int one_int = 1;
        const double alpha = 1.0, beta = 0.0;

        std::vector<double> rho_VU(this->LM->ParaV->nloc);

        for (int ik = 0; ik < kv.nks; ik++)
        {

            const int spin = kv.isk[ik];

            double* VU = new double[this->LM->ParaV->nloc];
            this->cal_VU_pot_mat_real(spin, false, VU);
            ModuleBase::timer::tick("DFTU", "cal_rho_VU");

#ifdef __MPI
            pdgemm_(&transT, &transN,
                    &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
                    &alpha, 
                    dm_gamma[spin].c, &one_int, &one_int, this->LM->ParaV->desc, 
                    VU, &one_int, &one_int, this->LM->ParaV->desc,
                    &beta,
                    &rho_VU[0], &one_int, &one_int, this->LM->ParaV->desc);
#endif

            delete[] VU;
            ModuleBase::timer::tick("DFTU", "cal_rho_VU");
            if (GlobalV::CAL_FORCE)  this->cal_force_gamma(&rho_VU[0], force_dftu);
            if (GlobalV::CAL_STRESS) this->cal_stress_gamma(&rho_VU[0], stress_dftu);
        } // ik
    }
    else
    {
        const char transN = 'N', transT = 'T';
        const int one_int = 1;
			const std::complex<double> alpha(1.0,0.0), beta(0.0,0.0);

        std::vector<std::complex<double>> rho_VU(this->LM->ParaV->nloc);

        for (int ik = 0; ik < kv.nks; ik++)
        {
            const int spin = kv.isk[ik];

            std::complex<double>* VU = new std::complex<double>[this->LM->ParaV->nloc];
            this->cal_VU_pot_mat_complex(spin, false, VU);
            ModuleBase::timer::tick("DFTU", "cal_rho_VU");

#ifdef __MPI
            pzgemm_(&transT, &transN,
                    &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
                    &alpha, 
                    dm_k[ik].c, &one_int, &one_int, this->LM->ParaV->desc, 
                    VU, &one_int, &one_int, this->LM->ParaV->desc,
                    &beta,
                    &rho_VU[0], &one_int, &one_int, this->LM->ParaV->desc);
#endif

            delete[] VU;
            ModuleBase::timer::tick("DFTU", "cal_rho_VU");

            if (GlobalV::CAL_FORCE)  cal_force_k (ik, &rho_VU[0], force_dftu, kv.kvec_d);
            if (GlobalV::CAL_STRESS) cal_stress_k(ik, &rho_VU[0], stress_dftu, kv.kvec_d);
        } // ik
    }

#ifdef __MPI
    if (GlobalV::CAL_FORCE)
    {
        Parallel_Reduce::reduce_double_pool(force_dftu.c, force_dftu.nr * force_dftu.nc);
    }
#endif

    if (GlobalV::CAL_STRESS)
    {
#ifdef __MPI
        Parallel_Reduce::reduce_double_pool(stress_dftu.c, stress_dftu.nr * stress_dftu.nc);
#endif

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (i > j) stress_dftu(i, j) = stress_dftu(j, i);
            }
        }

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                stress_dftu(i, j) *= GlobalC::ucell.lat0 / GlobalC::ucell.omega;
            }
        }
    }
    ModuleBase::timer::tick("DFTU", "force_stress");

    return;
}

void DFTU::cal_force_k(const int ik, 
                    const std::complex<double>* rho_VU, 
                    ModuleBase::matrix& force_dftu,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("DFTU", "cal_force_k");
    ModuleBase::timer::tick("DFTU", "cal_force_k");

    const char transN = 'N', transC = 'C';
    const int one_int = 1;
    const std::complex<double> zero(0.0,0.0), one(1.0,0.0);

    std::vector<std::complex<double>> dm_VU_dSm(this->LM->ParaV->nloc);
    std::vector<std::complex<double>> dSm_k(this->LM->ParaV->nloc);

    for (int dim = 0; dim < 3; dim++)
    {
        this->folding_matrix_k(ik, dim + 1, 0, &dSm_k[0], kvec_d);

#ifdef __MPI
        pzgemm_(&transN, &transC,
                &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
                &one, 
                &dSm_k[0], &one_int, &one_int, this->LM->ParaV->desc, 
                rho_VU, &one_int, &one_int, this->LM->ParaV->desc,
                &zero,
                &dm_VU_dSm[0], &one_int, &one_int, this->LM->ParaV->desc);
#endif

        for (int ir = 0; ir < this->LM->ParaV->nrow; ir++)
        {
            const int iwt1 = this->LM->ParaV->local2global_row(ir);
            const int iat1 = GlobalC::ucell.iwt2iat[iwt1];

            for (int ic = 0; ic < this->LM->ParaV->ncol; ic++)
            {
                const int iwt2 = this->LM->ParaV->local2global_col(ic);
                const int irc = ic * this->LM->ParaV->nrow + ir;

                if (iwt1 == iwt2) force_dftu(iat1, dim) += dm_VU_dSm[irc].real();

            } // end ic
        } // end ir

#ifdef __MPI
        pzgemm_(&transN, &transN,
                &GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
                &one, 
                &dSm_k[0], &one_int, &one_int, this->LM->ParaV->desc, 
                rho_VU, &one_int, &one_int, this->LM->ParaV->desc,
                &zero,
                &dm_VU_dSm[0], &one_int, &one_int, this->LM->ParaV->desc);
#endif

        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            const int NL = GlobalC::ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

            if (LC == -1) continue;
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                const int iat = GlobalC::ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it]) continue;
                    const int N = GlobalC::ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0) continue;

                        for (int m = 0; m < 2 * l + 1; m++)
                        {
                            for (int ipol = 0; ipol < GlobalV::NPOL; ipol++)
                            {
                                const int iwt = this->iatlnmipol2iwt[iat][l][n][m][ipol];
                                const int mu = this->LM->ParaV->global2local_row(iwt);
                                const int nu = this->LM->ParaV->global2local_col(iwt);
                                if (mu < 0 || nu < 0) continue;

                                force_dftu(iat, dim) += dm_VU_dSm[nu * this->LM->ParaV->nrow + mu].real();
                            }
                        } //
                    } // n
                } // l
            } // ia
        } // it
    } // end dim
    ModuleBase::timer::tick("DFTU", "cal_force_k");

    return;
}

void DFTU::cal_stress_k(const int ik,
                        const std::complex<double>* rho_VU,
                        ModuleBase::matrix& stress_dftu,
                        const std::vector<ModuleBase::Vector3<double>>& kvec_d)
{
    ModuleBase::TITLE("DFTU", "cal_stress_k");
    ModuleBase::timer::tick("DFTU", "cal_stress_k");
    const char transN = 'N';
    const int one_int = 1;
    const std::complex<double> minus_half(-0.5,0.0), zero(0.0,0.0), one(1.0,0.0);

    std::vector<std::complex<double>> dm_VU_sover(this->LM->ParaV->nloc);
    std::vector<std::complex<double>> dSR_k(this->LM->ParaV->nloc);

    for (int dim1 = 0; dim1 < 3; dim1++)
    {
        for (int dim2 = dim1; dim2 < 3; dim2++)
        {
            this->folding_matrix_k(ik, dim1 + 4, dim2, &dSR_k[0], kvec_d);

#ifdef __MPI
            pzgemm_(&transN, &transN,
					&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
					&minus_half, 
					rho_VU, &one_int, &one_int, this->LM->ParaV->desc, 
					&dSR_k[0], &one_int, &one_int, this->LM->ParaV->desc,
					&zero,
					&dm_VU_sover[0], &one_int, &one_int, this->LM->ParaV->desc);
#endif

            for (int ir = 0; ir < this->LM->ParaV->nrow; ir++)
            {
                const int iwt1 = this->LM->ParaV->local2global_row(ir);
                for (int ic = 0; ic < this->LM->ParaV->ncol; ic++)
                {
                    const int iwt2 = this->LM->ParaV->local2global_col(ic);
                    const int irc = ic * this->LM->ParaV->nrow + ir;

                    if (iwt1 == iwt2) stress_dftu(dim1, dim2) += 2.0 * dm_VU_sover[irc].real();
                } // end ic
            } // end ir

        } // end dim2
    } // end dim1
    ModuleBase::timer::tick("DFTU", "cal_stress_k");

    return;
}

void DFTU::cal_force_gamma(const double* rho_VU, ModuleBase::matrix& force_dftu)
{
    ModuleBase::TITLE("DFTU", "cal_force_gamma");
    ModuleBase::timer::tick("DFTU", "cal_force_gamma");
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const double one = 1.0, zero = 0.0, minus_one = -1.0;

    std::vector<double> dm_VU_dSm(this->LM->ParaV->nloc);

    for (int dim = 0; dim < 3; dim++)
    {
        double* tmp_ptr;
        if (dim == 0)      tmp_ptr = this->LM->DSloc_x;
        else if (dim == 1) tmp_ptr = this->LM->DSloc_y;
        else if (dim == 2) tmp_ptr = this->LM->DSloc_z;

#ifdef __MPI
        pdgemm_(&transN, &transT,
				&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
				&one, 
				tmp_ptr, &one_int, &one_int, this->LM->ParaV->desc, 
				rho_VU, &one_int, &one_int, this->LM->ParaV->desc,
				&zero,
				&dm_VU_dSm[0], &one_int, &one_int, this->LM->ParaV->desc);
#endif

        for (int ir = 0; ir < this->LM->ParaV->nrow; ir++)
        {
            const int iwt1 = this->LM->ParaV->local2global_row(ir);
            const int iat1 = GlobalC::ucell.iwt2iat[iwt1];

            for (int ic = 0; ic < this->LM->ParaV->ncol; ic++)
            {
                const int iwt2 = this->LM->ParaV->local2global_col(ic);
                const int irc = ic * this->LM->ParaV->nrow + ir;

                if (iwt1 == iwt2) force_dftu(iat1, dim) += dm_VU_dSm[irc];

            } // end ic
        } // end ir

#ifdef __MPI
        pdgemm_(&transN, &transT,
				&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
				&one, 
				tmp_ptr, &one_int, &one_int, this->LM->ParaV->desc, 
				rho_VU, &one_int, &one_int, this->LM->ParaV->desc,
				&zero,
				&dm_VU_dSm[0], &one_int, &one_int, this->LM->ParaV->desc);
#endif

        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            const int NL = GlobalC::ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

            if (LC == -1) continue;
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                const int iat = GlobalC::ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it]) continue;

                    const int N = GlobalC::ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0) continue;

                        // Calculate the local occupation number matrix
                        for (int m = 0; m < 2 * l + 1; m++)
                        {
                            for (int ipol = 0; ipol < GlobalV::NPOL; ipol++)
                            {
                                const int iwt = this->iatlnmipol2iwt[iat][l][n][m][ipol];
                                const int mu = this->LM->ParaV->global2local_row(iwt);
                                const int nu = this->LM->ParaV->global2local_col(iwt);
                                if (mu < 0 || nu < 0) continue;

                                force_dftu(iat, dim) += dm_VU_dSm[nu * this->LM->ParaV->nrow + mu];
                            }
                        } //
                    } // n
                } // l
            } // ia
        } // it

    } // end dim
    ModuleBase::timer::tick("DFTU", "cal_force_gamma");

    return;
}

void DFTU::cal_stress_gamma(const double* rho_VU, ModuleBase::matrix& stress_dftu)
{
    ModuleBase::TITLE("DFTU", "cal_stress_gamma");
    ModuleBase::timer::tick("DFTU", "cal_stress_gamma");
    const char transN = 'N';
    const int one_int = 1;
    const double zero = 0.0, minus_half = -0.5, one = 1.0;

    std::vector<double> dSR_gamma(this->LM->ParaV->nloc);
    std::vector<double> dm_VU_sover(this->LM->ParaV->nloc);

    for (int dim1 = 0; dim1 < 3; dim1++)
    {
        for (int dim2 = dim1; dim2 < 3; dim2++)
        {
            this->fold_dSR_gamma(dim1, dim2, &dSR_gamma[0]);

#ifdef __MPI
            pdgemm_(&transN, &transN,
					&GlobalV::NLOCAL, &GlobalV::NLOCAL, &GlobalV::NLOCAL,
					&minus_half, 
					rho_VU, &one_int, &one_int, this->LM->ParaV->desc, 
					&dSR_gamma[0], &one_int, &one_int, this->LM->ParaV->desc,
					&zero,
					&dm_VU_sover[0], &one_int, &one_int, this->LM->ParaV->desc);
#endif

            for (int ir = 0; ir < this->LM->ParaV->nrow; ir++)
            {
                const int iwt1 = this->LM->ParaV->local2global_row(ir);

                for (int ic = 0; ic < this->LM->ParaV->ncol; ic++)
                {
                    const int iwt2 = this->LM->ParaV->local2global_col(ic);
                    const int irc = ic * this->LM->ParaV->nrow + ir;

                    if (iwt1 == iwt2) stress_dftu(dim1, dim2) += 2.0 * dm_VU_sover[irc];
                } // end ic
            } // end ir

        } // end dim2
    } // end dim1
    ModuleBase::timer::tick("DFTU", "cal_stress_gamma");
    return;
}
} // namespace ModuleDFTU