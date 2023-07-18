#include "dftu.h"
#include "module_base/timer.h"
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
void DFTU::copy_locale()
{
    ModuleBase::TITLE("DFTU", "copy_locale");
    ModuleBase::timer::tick("DFTU", "copy_locale");

    for (int T = 0; T < GlobalC::ucell.ntype; T++)
    {
        if (orbital_corr[T] == -1)
            continue;

        for (int I = 0; I < GlobalC::ucell.atoms[T].na; I++)
        {
            const int iat = GlobalC::ucell.itia2iat(T, I);

            for (int l = 0; l < GlobalC::ucell.atoms[T].nwl + 1; l++)
            {
                const int N = GlobalC::ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (GlobalV::NSPIN == 4)
                    {
                        locale_save[iat][l][n][0] = locale[iat][l][n][0];
                    }
                    else if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
                    {
                        locale_save[iat][l][n][0] = locale[iat][l][n][0];
                        locale_save[iat][l][n][1] = locale[iat][l][n][1];
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DFTU", "copy_locale");
}

void DFTU::zero_locale()
{
    ModuleBase::TITLE("DFTU", "zero_locale");
    ModuleBase::timer::tick("DFTU", "zero_locale");

    for (int T = 0; T < GlobalC::ucell.ntype; T++)
    {
        if (orbital_corr[T] == -1) continue;

        for (int I = 0; I < GlobalC::ucell.atoms[T].na; I++)
        {
            const int iat = GlobalC::ucell.itia2iat(T, I);

            for (int l = 0; l < GlobalC::ucell.atoms[T].nwl + 1; l++)
            {
                const int N = GlobalC::ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (GlobalV::NSPIN == 4)
                    {
                        locale[iat][l][n][0].zero_out();
                    }
                    else if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
                    {
                        locale[iat][l][n][0].zero_out();
                        locale[iat][l][n][1].zero_out();
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DFTU", "zero_locale");
}

void DFTU::mix_locale(const double& mixing_beta)
{
    ModuleBase::TITLE("DFTU", "mix_locale");
    ModuleBase::timer::tick("DFTU", "mix_locale");

    double beta = mixing_beta;

    for (int T = 0; T < GlobalC::ucell.ntype; T++)
    {
        if (orbital_corr[T] == -1)
            continue;

        for (int I = 0; I < GlobalC::ucell.atoms[T].na; I++)
        {
            const int iat = GlobalC::ucell.itia2iat(T, I);

            for (int l = 0; l < GlobalC::ucell.atoms[T].nwl + 1; l++)
            {
                const int N = GlobalC::ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (GlobalV::NSPIN == 4)
                    {
                        locale[iat][l][n][0] = locale[iat][l][n][0]*beta + locale_save[iat][l][n][0]*(1.0-beta);
                    }
                    else if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 2)
                    {
                        locale[iat][l][n][0] = locale[iat][l][n][0] * beta + locale_save[iat][l][n][0] * (1.0-beta);
                        locale[iat][l][n][1] = locale[iat][l][n][1] * beta + locale_save[iat][l][n][1] * (1.0-beta);
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DFTU", "mix_locale");
}

void DFTU::cal_occup_m_k(const int iter, 
                        std::vector<ModuleBase::ComplexMatrix> &dm_k,
                        const K_Vectors& kv,
                        const double& mixing_beta)
{
    ModuleBase::TITLE("DFTU", "cal_occup_m_k");
    ModuleBase::timer::tick("DFTU", "cal_occup_m_k");

    this->copy_locale();
    this->zero_locale();

    //=================Part 1======================
    // call SCALAPACK routine to calculate the product of the S and density matrix
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const std::complex<double> beta(0.0,0.0), alpha(1.0,0.0);

    std::vector<std::complex<double>> srho(this->LM->ParaV->nloc);
    std::vector<std::complex<double>> Sk(this->LM->ParaV->nloc);

    for (int ik = 0; ik < kv.nks; ik++)
    {
        // srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_k(iw,nu)
        this->folding_matrix_k(ik, 0, 0, &Sk[0], kv.kvec_d);

#ifdef __MPI
        pzgemm_(&transN,
                &transT,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &alpha,
                &Sk[0],
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                dm_k[ik].c,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                &beta,
                &srho[0],
                &one_int,
                &one_int,
                this->LM->ParaV->desc);
#endif

        const int spin = kv.isk[ik];
        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            const int NL = GlobalC::ucell.atoms[it].nwl + 1;
            const int LC = orbital_corr[it];

            if (LC == -1)
                continue;

            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                const int iat = GlobalC::ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it])
                        continue;

                    const int N = GlobalC::ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        // if(!Yukawa && n!=0) continue;
                        if (n != 0)
                            continue;

                        // Calculate the local occupation number matrix
                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < GlobalV::NPOL; ipol0++)
                            {
                                const int iwt0 = this->iatlnmipol2iwt[iat][l][n][m0][ipol0];
                                const int mu = this->LM->ParaV->global2local_row(iwt0);
                                const int mu_prime = this->LM->ParaV->global2local_col(iwt0);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < GlobalV::NPOL; ipol1++)
                                    {
                                        const int iwt1 = this->iatlnmipol2iwt[iat][l][n][m1][ipol1];
                                        const int nu = this->LM->ParaV->global2local_col(iwt1);
                                        const int nu_prime = this->LM->ParaV->global2local_row(iwt1);

                                        const int irc = nu * this->LM->ParaV->nrow + mu;
                                        const int irc_prime = mu_prime * this->LM->ParaV->nrow + nu_prime;

                                        const int m0_all = m0 + ipol0 * (2 * l + 1);
                                        const int m1_all = m1 + ipol1 * (2 * l + 1);

                                        if ((nu >= 0) && (mu >= 0))
                                            locale[iat][l][n][spin](m0_all, m1_all) += (srho[irc]).real() / 4.0;

                                        if ((nu_prime >= 0) && (mu_prime >= 0))
                                            locale[iat][l][n][spin](m0_all, m1_all)
                                                += (std::conj(srho[irc_prime])).real() / 4.0;
                                    } // ipol1
                                } // m1
                            } // ipol0
                        } // m0
                    } // end n
                } // end l
            } // end ia
        } // end it
    } // ik

    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        const int NL = GlobalC::ucell.atoms[it].nwl + 1;
        const int LC = orbital_corr[it];

        if (LC == -1)
            continue;

        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            const int iat = GlobalC::ucell.itia2iat(it, ia);

            for (int l = 0; l < NL; l++)
            {
                if (l != orbital_corr[it])
                    continue;

                const int N = GlobalC::ucell.atoms[it].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    // if(!Yukawa && n!=0) continue;
                    if (n != 0)
                        continue;
                        // set the local occupation mumber matrix of spin up and down zeros

#ifdef __MPI
                    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
                    {
                        ModuleBase::matrix temp(locale[iat][l][n][0]);
                        MPI_Allreduce(&temp(0, 0),
                                      &locale[iat][l][n][0](0, 0),
                                      (2 * l + 1) * GlobalV::NPOL * (2 * l + 1) * GlobalV::NPOL,
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
                    }
                    else if (GlobalV::NSPIN == 2)
                    {
                        ModuleBase::matrix temp0(locale[iat][l][n][0]);
                        MPI_Allreduce(&temp0(0, 0),
                                      &locale[iat][l][n][0](0, 0),
                                      (2 * l + 1) * (2 * l + 1),
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);

                        ModuleBase::matrix temp1(locale[iat][l][n][1]);
                        MPI_Allreduce(&temp1(0, 0),
                                      &locale[iat][l][n][1](0, 0),
                                      (2 * l + 1) * (2 * l + 1),
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
                    }
#endif

                    // for the case spin independent calculation
                    switch (GlobalV::NSPIN)
                    {
                    case 1:
                        locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
                        locale[iat][l][n][0] *= 0.5;
                        locale[iat][l][n][1] += locale[iat][l][n][0];
                        break;

                    case 2:
                        for (int is = 0; is < GlobalV::NSPIN; is++)
                            locale[iat][l][n][is] += transpose(locale[iat][l][n][is]);
                        break;

                    case 4: // SOC
                        locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
                        break;

                    default:
                        std::cout << "Not supported NSPIN parameter" << std::endl;
                        exit(0);
                    }
                } // end n
            } // end l
        } // end ia
    } // end it

    if(mixing_dftu && initialed_locale)
    {
        this->mix_locale(mixing_beta);
    }

    this->initialed_locale = true;
    ModuleBase::timer::tick("DFTU", "cal_occup_m_k");
    return;
}

void DFTU::cal_occup_m_gamma(const int iter, std::vector<ModuleBase::matrix> &dm_gamma, const double& mixing_beta)
{
    ModuleBase::TITLE("DFTU", "cal_occup_m_gamma");
    ModuleBase::timer::tick("DFTU", "cal_occup_m_gamma");
    this->copy_locale();
    this->zero_locale();

    //=================Part 1======================
    // call PBLAS routine to calculate the product of the S and density matrix
    const char transN = 'N', transT = 'T';
    const int one_int = 1;
    const double alpha = 1.0, beta = 0.0;

    std::vector<double> srho(this->LM->ParaV->nloc);
    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        // srho(mu,nu) = \sum_{iw} S(mu,iw)*dm_gamma(iw,nu)

#ifdef __MPI
        pdgemm_(&transN,
                &transT,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &GlobalV::NLOCAL,
                &alpha,
                this->LM->Sloc.data(),
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                dm_gamma[is].c,
                &one_int,
                &one_int,
                this->LM->ParaV->desc,
                &beta,
                &srho[0],
                &one_int,
                &one_int,
                this->LM->ParaV->desc);
#endif

        for (int it = 0; it < GlobalC::ucell.ntype; it++)
        {
            const int NL = GlobalC::ucell.atoms[it].nwl + 1;
            if (orbital_corr[it] == -1)
                continue;
            for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
            {
                const int iat = GlobalC::ucell.itia2iat(it, ia);

                for (int l = 0; l < NL; l++)
                {
                    if (l != orbital_corr[it])
                        continue;

                    const int N = GlobalC::ucell.atoms[it].l_nchi[l];

                    for (int n = 0; n < N; n++)
                    {
                        if (n != 0)
                            continue;

                        // Calculate the local occupation number matrix
                        for (int m0 = 0; m0 < 2 * l + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < GlobalV::NPOL; ipol0++)
                            {
                                const int iwt0 = this->iatlnmipol2iwt[iat][l][n][m0][ipol0];
                                const int mu = this->LM->ParaV->global2local_row(iwt0);
                                const int mu_prime = this->LM->ParaV->global2local_col(iwt0);

                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < GlobalV::NPOL; ipol1++)
                                    {
                                        const int iwt1 = this->iatlnmipol2iwt[iat][l][n][m1][ipol1];
                                        const int nu = this->LM->ParaV->global2local_col(iwt1);
                                        const int nu_prime = this->LM->ParaV->global2local_row(iwt1);

                                        const int irc = nu * this->LM->ParaV->nrow + mu;
                                        const int irc_prime = mu_prime * this->LM->ParaV->nrow + nu_prime;

                                        if ((nu >= 0) && (mu >= 0))
                                        {
                                            int m0_all = m0 + (2 * l + 1) * ipol0;
                                            int m1_all = m0 + (2 * l + 1) * ipol1;

                                            locale[iat][l][n][is](m0, m1) += srho[irc] / 4.0;
                                        }

                                        if ((nu_prime >= 0) && (mu_prime >= 0))
                                        {
                                            int m0_all = m0 + (2 * l + 1) * ipol0;
                                            int m1_all = m0 + (2 * l + 1) * ipol1;

                                            locale[iat][l][n][is](m0, m1) += srho[irc_prime] / 4.0;
                                        }
                                    }
                                }
                            }
                        }

                        ModuleBase::matrix temp(locale[iat][l][n][is]);

#ifdef __MPI
                        MPI_Allreduce(&temp(0, 0),
                                      &locale[iat][l][n][is](0, 0),
                                      (2 * l + 1) * GlobalV::NPOL * (2 * l + 1) * GlobalV::NPOL,
                                      MPI_DOUBLE,
                                      MPI_SUM,
                                      MPI_COMM_WORLD);
#endif

                        // for the case spin independent calculation
                        switch (GlobalV::NSPIN)
                        {
                        case 1:
                            locale[iat][l][n][0] += transpose(locale[iat][l][n][0]);
                            locale[iat][l][n][0] *= 0.5;
                            locale[iat][l][n][1] += locale[iat][l][n][0];
                            break;

                        case 2:
                            locale[iat][l][n][is] += transpose(locale[iat][l][n][is]);
                            break;

                        default:
                            std::cout << "Not supported NSPIN parameter" << std::endl;
                            exit(0);
                        }

                    } // end for(n)
                } // L
            } // ia
        } // it
    } // is

    if(mixing_dftu && initialed_locale)
    {
        this->mix_locale(mixing_beta);
    }

    this->initialed_locale = true;
    ModuleBase::timer::tick("DFTU", "cal_occup_m_gamma");
    return;
}
} // namespace ModuleDFTU