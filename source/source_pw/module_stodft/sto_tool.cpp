#include "sto_tool.h"

#include "source_base/math_chebyshev.h"
#include "source_base/parallel_device.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include <vector>

template <typename FPTYPE, typename Device>
void check_che_op<FPTYPE, Device>::operator()(const int& nche_in,
                                              const double& try_emin,
                                              const double& try_emax,
                                              const int& nbands_sto,
                                              K_Vectors* p_kv,
                                              Stochastic_WF<std::complex<FPTYPE>, Device>* p_stowf,
                                              hamilt::HamiltSdftPW<std::complex<FPTYPE>, Device>* p_hamilt_sto)
{
    //------------------------------
    //      Convergence test
    //------------------------------
    bool change = false;
    const int nk = p_kv->get_nks();
    ModuleBase::Chebyshev<FPTYPE, Device> chetest(nche_in);
    int ntest0 = 5;
    *p_hamilt_sto->emax = try_emax;
    *p_hamilt_sto->emin = try_emin;
    // if (PARAM.inp.nbands > 0)
    // {
    //     double tmpemin = 1e10;
    //     for (int ik = 0; ik < nk; ++ik)
    //     {
    //         tmpemin = std::min(tmpemin, this->pelec->ekb(ik, PARAM.inp.nbands - 1));
    //     }
    //     *p_hamilt_sto->emin = tmpemin;
    // }
    // else
    // {
    //     *p_hamilt_sto->emin = 0;
    // }
    for (int ik = 0; ik < nk; ++ik)
    {
        p_hamilt_sto->updateHk(ik);
        const int npw = p_kv->ngk[ik];
        std::complex<FPTYPE>* pchi = nullptr;
        psi::Psi<std::complex<FPTYPE>, Device> randchi_d;
        if (nbands_sto == 0) // For case: PARAM.inp.nbands_sto = "all"
        {
            randchi_d.resize(1, 1, npw);
        }
        int ntest = std::min(ntest0, p_stowf->nchip[ik]);
        for (int i = 0; i < ntest; ++i)
        {
            if (nbands_sto == 0)
            {
                std::vector<std::complex<FPTYPE>> randchi(npw);   
                for (int ig = 0; ig < npw; ++ig)
                {
                    FPTYPE rr = std::rand() / FPTYPE(RAND_MAX);
                    FPTYPE arg = std::rand() / FPTYPE(RAND_MAX);
                    randchi[ig] = std::complex<FPTYPE>(rr * cos(arg), rr * sin(arg));
                }
                syncmem_complex_h2d_op()(randchi_d.get_pointer(), randchi.data(), npw);
                pchi = randchi_d.get_pointer();
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
                auto hchi_norm = std::bind(&hamilt::HamiltSdftPW<std::complex<FPTYPE>, Device>::hPsi_norm,
                                           p_hamilt_sto,
                                           std::placeholders::_1,
                                           std::placeholders::_2,
                                           std::placeholders::_3);
                converge = chetest.checkconverge(hchi_norm,
                                                 pchi,
                                                 npw,
                                                 p_stowf->npwx,
                                                 *p_hamilt_sto->emax,
                                                 *p_hamilt_sto->emin,
                                                 2.0);

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
            MPI_Allreduce(MPI_IN_PLACE, p_hamilt_sto->emax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, p_hamilt_sto->emin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
            GlobalV::ofs_running << "New Emax " << *p_hamilt_sto->emax << " Ry; new Emin " << *p_hamilt_sto->emin
                                 << " Ry" << std::endl;
            change = false;
        }
    }
}

template <typename FPTYPE, typename Device>
psi::Psi<std::complex<FPTYPE>, Device>* gatherchi_op<FPTYPE, Device>::operator()(
    psi::Psi<std::complex<FPTYPE>, Device>& chi,
    psi::Psi<std::complex<FPTYPE>, Device>& chi_all,
    const int& npwx,
    int* nrecv_sto,
    int* displs_sto,
    const int perbands_sto)
{
    psi::Psi<std::complex<FPTYPE>, Device>* p_chi;
    p_chi = &chi;
#ifdef __MPI
    if (PARAM.inp.bndpar > 1)
    {
        p_chi = &chi_all;
        ModuleBase::timer::tick("sKG", "bands_gather");
        Parallel_Common::gatherv_dev<std::complex<FPTYPE>, Device>(chi.get_pointer(),
                                                                   perbands_sto * npwx,
                                                                   chi_all.get_pointer(),
                                                                   nrecv_sto,
                                                                   displs_sto,
                                                                   BP_WORLD);
        ModuleBase::timer::tick("sKG", "bands_gather");
    }
#endif
    return p_chi;
}

template struct check_che_op<double, base_device::DEVICE_CPU>;
#ifdef __ENABLE_FLOAT_FFTW
template struct check_che_op<float, base_device::DEVICE_CPU>;
#endif
template struct gatherchi_op<double, base_device::DEVICE_CPU>;
template struct gatherchi_op<float, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template struct check_che_op<double, base_device::DEVICE_GPU>;
#ifdef __ENABLE_FLOAT_FFTW
template struct check_che_op<float, base_device::DEVICE_GPU>;
#endif
template struct gatherchi_op<double, base_device::DEVICE_GPU>;
template struct gatherchi_op<float, base_device::DEVICE_GPU>;
#endif