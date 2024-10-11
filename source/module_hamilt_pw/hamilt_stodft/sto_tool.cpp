#include "sto_tool.h"

#include "module_parameter/parameter.h"
#include "module_base/math_chebyshev.h"
#include "module_base/timer.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include <vector>

void check_che(const int& nche_in,
               const double& try_emin,
               const double& try_emax,
               const int& nbands_sto,
               K_Vectors* p_kv,
               Stochastic_WF* p_stowf,
               hamilt::Hamilt<std::complex<double>>* p_hamilt,
               Stochastic_hchi& stohchi)
{
    //------------------------------
    //      Convergence test
    //------------------------------
    bool change = false;
    const int nk = p_kv->get_nks();
    ModuleBase::Chebyshev<double> chetest(nche_in);
    int ntest0 = 5;
    *stohchi.Emax = try_emax;
    *stohchi.Emin = try_emin;
    // if (PARAM.inp.nbands > 0)
    // {
    //     double tmpemin = 1e10;
    //     for (int ik = 0; ik < nk; ++ik)
    //     {
    //         tmpemin = std::min(tmpemin, this->pelec->ekb(ik, PARAM.inp.nbands - 1));
    //     }
    //     stohchi.Emin = tmpemin;
    // }
    // else
    // {
    //     stohchi.Emin = 0;
    // }
    for (int ik = 0; ik < nk; ++ik)
    {
        p_hamilt->updateHk(ik);
        stohchi.current_ik = ik;
        const int npw = p_kv->ngk[ik];
        std::complex<double>* pchi = nullptr;
        std::vector<std::complex<double>> randchi;
        int ntest = std::min(ntest0, p_stowf->nchip[ik]);
        for (int i = 0; i < ntest; ++i)
        {
            if (nbands_sto == 0)
            {
                randchi.resize(npw);
                pchi = &randchi[0];
                for (int ig = 0; ig < npw; ++ig)
                {
                    double rr = std::rand() / double(RAND_MAX);
                    double arg = std::rand() / double(RAND_MAX);
                    pchi[ig] = std::complex<double>(rr * cos(arg), rr * sin(arg));
                }
            }
            else if (PARAM.inp.nbands > 0)
            {
                pchi = &p_stowf->chiortho[0](ik, i, 0);
            }
            else
            {
                pchi = &p_stowf->chi0[0](ik, i, 0);
            }
            while (true)
            {
                bool converge;
                auto hchi_norm = std::bind(&Stochastic_hchi::hchi_norm,
                                           &stohchi,
                                           std::placeholders::_1,
                                           std::placeholders::_2,
                                           std::placeholders::_3);
                converge = chetest.checkconverge(hchi_norm, pchi, npw, *stohchi.Emax, *stohchi.Emin, 2.0);

                if (!converge)
                {
                    change = true;
                }
                else
                {
                    break;
                }
            }
        }

        if (ik == nk - 1)
        {
#ifdef __MPI
            MPI_Allreduce(MPI_IN_PLACE, stohchi.Emax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, stohchi.Emin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
            GlobalV::ofs_running << "New Emax " << *stohchi.Emax << " Ry; new Emin " << *stohchi.Emin << " Ry"
                                 << std::endl;
            change = false;
        }
    }
}

void convert_psi(const psi::Psi<std::complex<double>>& psi_in, psi::Psi<std::complex<float>>& psi_out)
{
    psi_in.fix_k(0);
    psi_out.fix_k(0);
    for (int i = 0; i < psi_in.size(); ++i)
    {
        psi_out.get_pointer()[i] = static_cast<std::complex<float>>(psi_in.get_pointer()[i]);
    }
    return;
}

psi::Psi<std::complex<float>>* gatherchi(psi::Psi<std::complex<float>>& chi,
                                         psi::Psi<std::complex<float>>& chi_all,
                                         const int& npwx,
                                         int* nrecv_sto,
                                         int* displs_sto,
                                         const int perbands_sto)
{
    psi::Psi<std::complex<float>>* p_chi;
    p_chi = &chi;
#ifdef __MPI
    if (PARAM.inp.bndpar > 1)
    {
        p_chi = &chi_all;
        ModuleBase::timer::tick("sKG", "bands_gather");
        MPI_Allgatherv(chi.get_pointer(),
                       perbands_sto * npwx,
                       MPI_COMPLEX,
                       chi_all.get_pointer(),
                       nrecv_sto,
                       displs_sto,
                       MPI_COMPLEX,
                       PARAPW_WORLD);
        ModuleBase::timer::tick("sKG", "bands_gather");
    }
#endif
    return p_chi;
}