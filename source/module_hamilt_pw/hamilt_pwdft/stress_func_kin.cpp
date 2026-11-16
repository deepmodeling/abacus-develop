#include "stress_func.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/fs_kin_tools.h"

#if defined(__CUDA)
#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>

template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_kin_cuda(ModuleBase::matrix& sigma, const ModuleBase::matrix& wg,
                                                ModuleSymmetry::Symmetry* symmetry, K_Vectors* kpoints,
                                                ModulePW::PW_Basis_K* basis, const UnitCell& ucell_in,
                                                const psi::Psi<complex<FPTYPE>, Device>* wavefunctions)
{
    static_assert(std::is_same<FPTYPE, double>::value, "Kinetic stress uses double precision");
    ModuleBase::timer::tick("Stress_Func", "stress_kin");
    const int npol = ucell_in.get_npol();
    const int npwx = wavefunctions->get_nbasis() / npol;
    const double tpiba = ModuleBase::TWO_PI / ucell_in.lat0;
    const int max_blocks = hamilt::kinetic_stress_blocks_gpu(basis->npwk_max);
    auto weights = ct::TensorMap(wg.c, ct::DataType::DT_DOUBLE, ct::DeviceType::CpuDevice, {wg.nr, wg.nc})
                       .to_device<ct::DEVICE_GPU>();
    ct::Tensor partial(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, {max_blocks, 6});
    std::vector<double> host_partial(6 * max_blocks);
    double sums[6] = {};
    for (int ik = 0; ik < kpoints->get_nks(); ++ik)
    {
        const int npw = kpoints->ngk[ik];
        if (npw == 0) { continue; }
        const int blocks = hamilt::kinetic_stress_blocks_gpu(npw);
        hamilt::kinetic_stress_gpu(npw, npwx, PARAM.inp.nbands, npol, tpiba,
                                   basis->erf_ecut, basis->erf_height, basis->erf_sigma,
                                   basis->template get_gcar_data<double>() + 3 * ik * basis->npwk_max,
                                   basis->template get_kvec_c_data<double>() + 3 * ik,
                                   weights.data<double>() + ik * wg.nc,
                                   &wavefunctions[0](ik, 0, 0), partial.data<double>());
        base_device::memory::synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(
            nullptr, nullptr, host_partial.data(), partial.data<double>(), 6 * blocks);
        for (int block = 0; block < blocks; ++block)
        {
            for (int component = 0; component < 6; ++component)
            {
                sums[component] += host_partial[6 * block + component];
            }
        }
    }
    const double scale = (PARAM.inp.gamma_only ? 2.0 : 1.0) * ModuleBase::e2 / ucell_in.omega;
    int component = 0;
    for (int row = 0; row < 3; ++row)
    {
        for (int column = 0; column <= row; ++column)
        {
            sigma(row, column) = sigma(column, row) = scale * sums[component++];
        }
    }
    Parallel_Reduce::reduce_all(sigma.c, 9);
    if (ModuleSymmetry::Symmetry::symm_flag == 1)
    {
        symmetry->symmetrize_mat3(sigma, ucell_in.lat);
    }
    ModuleBase::timer::tick("Stress_Func", "stress_kin");
}
#endif

//calculate the kinetic stress in PW base
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_kin(ModuleBase::matrix& sigma,
                                             const ModuleBase::matrix& wg,
                                             ModuleSymmetry::Symmetry* p_symm,
                                             K_Vectors* p_kv,
                                             ModulePW::PW_Basis_K* wfc_basis,
                                             const UnitCell& ucell_in,
                                             const psi::Psi<complex<FPTYPE>, Device>* psi_in)
{
    ModuleBase::TITLE("Stress_Func","stress_kin");
	ModuleBase::timer::tick("Stress_Func","stress_kin");
	this->ucell = &ucell_in;
	hamilt::FS_Kin_tools<FPTYPE, Device> kin_tool(*this->ucell, p_kv, wfc_basis, wg);
    for (int ik = 0; ik < wfc_basis->nks; ++ik)
    {
        int nbands_occ = wg.nc;
        while (wg(ik, nbands_occ - 1) == 0.0)
        {
            nbands_occ--;
            if (nbands_occ == 0)
            {
                break;
            }
        }
        const int npm = nbands_occ;
        kin_tool.cal_gk(ik);
        kin_tool.cal_stress_kin(ik, npm, true, &psi_in[0](ik, 0, 0));
    }
    kin_tool.symmetrize_stress(p_symm, sigma);
		
	ModuleBase::timer::tick("Stress_Func","stress_kin");
	return;
}

template class Stress_Func<double, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, base_device::DEVICE_GPU>;
#endif
