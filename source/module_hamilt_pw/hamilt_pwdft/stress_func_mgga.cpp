#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "stress_func.h"

// calculate the Pulay term of mGGA stress correction in PW
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_mgga(ModuleBase::matrix& sigma,
                                              const ModuleBase::matrix& wg,
                                              const ModuleBase::matrix& v_ofk,
                                              const Charge* const chr,
                                              K_Vectors* p_kv,
                                              ModulePW::PW_Basis_K* wfc_basis,
                                              const psi::Psi<complex<FPTYPE>>* psi_in)
{
    ModuleBase::timer::tick("Stress_Func", "stress_mgga");

    if (GlobalV::NSPIN == 4)
        ModuleBase::WARNING_QUIT("stress_mgga", "noncollinear stress + mGGA not implemented");

    int current_spin = 0;
    const int nrxx = wfc_basis->nrxx;

    std::complex<FPTYPE>* psi;
    int ipol2xy[3][3];
    FPTYPE sigma_mgga[3][3];

    std::complex<FPTYPE>* gradwfc = new std::complex<FPTYPE>[nrxx * 3];
    ModuleBase::GlobalFunc::ZEROS(gradwfc, nrxx * 3);

    FPTYPE** crosstaus = new FPTYPE*[GlobalV::NSPIN];
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        crosstaus[is] = new FPTYPE[nrxx * 6];
        ModuleBase::GlobalFunc::ZEROS(crosstaus[is], GlobalV::NSPIN);
    }

    for (int ik = 0; ik < p_kv->nks; ik++)
    {
        if (GlobalV::NSPIN == 2)
            current_spin = p_kv->isk[ik];
        const int npw = p_kv->ngk[ik];
        psi = new complex<FPTYPE>[npw];

        for (int ibnd = 0; ibnd < GlobalV::NBANDS; ibnd++)
        {
            const FPTYPE w1 = wg(ik, ibnd) / GlobalC::ucell.omega;
            const std::complex<FPTYPE>* ppsi = nullptr;
            ppsi = &(psi_in[0](ik, ibnd, 0));
            for (int ig = 0; ig < npw; ig++)
            {
                psi[ig] = ppsi[ig];
            }
            XC_Functional::grad_wfc(psi, ik, gradwfc, wfc_basis, GlobalC::ucell.tpiba);

            int ipol = 0;
            for (int ix = 0; ix < 3; ix++)
            {
                for (int iy = 0; iy < ix + 1; iy++)
                {
                    ipol2xy[ix][iy] = ipol;
                    ipol2xy[iy][ix] = ipol;
                    for (int ir = 0; ir < nrxx; ir++)
                    {
                        crosstaus[current_spin][ipol*nrxx + ir] += 2.0 * w1
                                                             * (gradwfc[ix*nrxx + ir].real() * gradwfc[iy*nrxx + ir].real()
                                                                + gradwfc[ix*nrxx + ir].imag() * gradwfc[iy*nrxx + ir].imag());
                    }
                    ipol += 1;
                }
            }
        } // band loop
        delete[] psi;
    } // k loop

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ipol = 0; ipol < 6; ++ipol)
        {
            chr->reduce_diff_pools(crosstaus[is] + ipol * nrxx);
        }
    }

    delete[] gradwfc;

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        for (int ix = 0; ix < 3; ix++)
        {
            for (int iy = 0; iy < 3; iy++)
            {
                FPTYPE delta = 0.0;
                if (ix == iy)
                    delta = 1.0;
                sigma_mgga[ix][iy] = 0.0;
                for (int ir = 0; ir < nrxx; ir++)
                {
                    FPTYPE x
                        = v_ofk(is, ir) * (chr->kin_r[is][ir] * delta + crosstaus[is][ipol2xy[ix][iy] * nrxx + ir]);
                    sigma_mgga[ix][iy] += x;
                }
            }
        }
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        delete[] crosstaus[is];
    }
    delete[] crosstaus;

#ifdef __MPI
    for (int l = 0; l < 3; l++)
    {
        for (int m = 0; m < 3; m++)
        {
            Parallel_Reduce::reduce_double_pool(sigma_mgga[l][m]);
        }
    }
#endif
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            sigma(i, j) += sigma_mgga[i][j] / wfc_basis->nxyz;
        }
    }
    ModuleBase::timer::tick("Stress_Func", "stress_mgga");
    return;
}

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif