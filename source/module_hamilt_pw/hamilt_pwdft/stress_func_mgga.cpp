#include "module_base/timer.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "stress_func.h"

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

// calculate the Pulay term of mGGA stress correction in PW
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_mgga(ModuleBase::matrix& sigma,
                                              const ModuleBase::matrix& wg,
                                              const ModuleBase::matrix& v_ofk,
                                              const Charge* const chr,
                                              K_Vectors* p_kv,
                                              ModulePW::PW_Basis_K* wfc_basis,
                                              const psi::Psi<complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::timer::tick("Stress_Func", "stress_mgga");

    if (GlobalV::NSPIN == 4)
        ModuleBase::WARNING_QUIT("stress_mgga", "noncollinear stress + mGGA not implemented");

    int current_spin = 0;
    const int nrxx = wfc_basis->nrxx;

    std::complex<FPTYPE>* psi;
    int ipol2xy[3][3] = {{0, 1, 3}, {1, 2, 4}, {3, 4, 5}};
    FPTYPE sigma_mgga[3][3];

    using ct_Device = typename ct::PsiToContainer<Device>::type;

    auto gradwfc = ct::Tensor(
        ct::DataTypeToEnum<std::complex<FPTYPE>>::value,
        ct::DeviceTypeToEnum<ct_Device>::value,
        {nrxx * 3});
    auto crosstaus = ct::Tensor(
        ct::DataTypeToEnum<FPTYPE>::value,
        ct::DeviceTypeToEnum<ct_Device>::value,
        {GlobalV::NSPIN, nrxx * 6});
    crosstaus.zero(); // Must be zeroed out 

    auto cal_stress_mgga_solver = hamilt::cal_stress_mgga_op<std::complex<FPTYPE>, Device>();
    for (int ik = 0; ik < p_kv->nks; ik++)
    {
        if (GlobalV::NSPIN == 2)
            current_spin = p_kv->isk[ik];
        const int npw = p_kv->ngk[ik];

        for (int ibnd = 0; ibnd < GlobalV::NBANDS; ibnd++)
        {
            const FPTYPE w1 = wg(ik, ibnd) / GlobalC::ucell.omega;
            const std::complex<FPTYPE>* psi = &psi_in[0](ik, ibnd, 0);
            XC_Functional::grad_wfc<std::complex<FPTYPE>, Device>(ik, GlobalC::ucell.tpiba, wfc_basis, psi, gradwfc.data<std::complex<FPTYPE>>());
            cal_stress_mgga_solver(
                current_spin, nrxx, w1, gradwfc.data<std::complex<FPTYPE>>(), crosstaus.data<FPTYPE>());
        } // band loop
        // delete[] psi;
    } // k loop
    auto crosstaus_host = crosstaus.to_device<ct::DEVICE_CPU>();
    auto crosstaus_pack = crosstaus_host.accessor<FPTYPE, 2>();
#ifdef __MPI
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        for (int ipol = 0; ipol < 6; ++ipol)
        {
            chr->reduce_diff_pools(&crosstaus_pack[is][ipol * nrxx]);
        }
    }
#endif

    for(int ix = 0; ix < 3; ix++)
    {
        for(int iy = 0; iy < 3; iy++)
        {
            sigma_mgga[ix][iy] = 0.0;
        }
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        for (int ix = 0; ix < 3; ix++)
        {
            for (int iy = 0; iy < 3; iy++)
            {
                FPTYPE delta = 0.0;
                if (ix == iy)
                    delta = 1.0;
                for (int ir = 0; ir < nrxx; ir++)
                {
                    FPTYPE x
                        = v_ofk(is, ir) * (chr->kin_r[is][ir] * delta + crosstaus_pack[is][ipol2xy[ix][iy] * nrxx + ir]);
                    sigma_mgga[ix][iy] += x;
                }
            }
        }
    }

#ifdef __MPI
    for (int l = 0; l < 3; l++)
    {
        for (int m = 0; m < 3; m++)
        {
            Parallel_Reduce::reduce_pool(sigma_mgga[l][m]);
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