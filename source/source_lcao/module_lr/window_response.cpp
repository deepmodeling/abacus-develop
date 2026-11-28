#include "window_response.h"
#include "window_response_math.h"
#include "gradient_inputs.h"
#include "hamilt_casida.h"
#include "hamilt_ulr.h"
#include "utils/lr_util.h"
#include "source_base/timer.h"
#include <algorithm>
#include <complex>
#include <stdexcept>

namespace LR
{
namespace
{
std::vector<container::Tensor> transform_window_matrix(const GradientInputs<double>& inputs,
    const psi::Psi<double>& c, const Parallel_Blkcyclic& pm, const double* matrix,
    const LR_Util::MO_TYPE type, const int spin)
{
#ifdef __MPI
    auto result = cal_dm_trans_pblas(matrix, pm, c, inputs.pc, inputs.nbasis,
        inputs.nocc[spin], inputs.nvirt[spin], inputs.pmat.para_2d(), 1.0, type);
    for (auto& tensor : result)
    {
        LR_Util::matsym(tensor.data<double>(), inputs.nbasis, inputs.pmat.para_2d());
    }
#else
    auto result = cal_dm_trans_blas(matrix, c, inputs.nocc[spin], inputs.nvirt[spin], 1.0, type);
    for (auto& tensor : result)
    {
        LR_Util::matsym(tensor.data<double>(), inputs.nbasis);
    }
#endif
    return result;
}
}

template <>
WindowResponse<double> build_window_response(const GradientInputs<double>& inputs,
    const container::Tensor& X, const int nstates, const int nactive, const int ispin,
    const std::string& in_dir)
{
    ModuleBase::timer::start("LR", "build_window_response");
    WindowResponse<double> result;
    const int nv = inputs.nvirt[0];
    const int no = inputs.nocc[0];
    if (nactive == nv)
    {
        ModuleBase::timer::end("LR", "build_window_response");
        return result;
    }
    if (inputs.nk != 1)
    {
        throw std::runtime_error("Window response currently requires real Gamma orbitals");
    }
    const auto c = LR_Util::get_psi_spin(inputs.psi_ks, 0, inputs.nk);
    const Parallel_Blkcyclic& px = inputs.px[0];
    Parallel_Blkcyclic pvv;
    LR_Util::setup_2d_division(pvv, px.nb(), nv, nv
#ifdef __MPI
        , px.ictxt()
#endif
    );
    std::weak_ptr<PotLRBase> pot_lr = inputs.pot[ispin];
    HamiltLR<double> full_casida(inputs.xc_kernel, inputs.nspin, inputs.nbasis,
        inputs.nocc, inputs.nvirt, inputs.ucell, inputs.orb_cutoff, inputs.gd,
        inputs.psi_ks, inputs.eig_ks,
#ifdef __EXX
        inputs.exx_lri, inputs.hybrid_alpha,
#endif
        pot_lr, inputs.kv, inputs.px, inputs.pc, inputs.pmat,
        inputs.spin_types[ispin], in_dir, inputs.out_dir);
    const int batch_size = nstates * inputs.nloc;
    const int batch_storage = std::max(batch_size, 1);
    std::vector<double> AX(batch_storage, 0.0);
    const double empty_block = 0.0;
    const double* const x_data = X.data<double>() == nullptr ? &empty_block : X.data<double>();
    full_casida.hPsi(x_data, AX.data(), inputs.nloc, nstates);

    const int pairs = nv * no;
    const int vv_storage = std::max<int>(pvv.size(), 1);
    const double* const virtual_energies = inputs.eig_ks.c + no;
    for (int state = 0; state < nstates; ++state)
    {
        const int offset = state * inputs.nloc;
        const double* x = x_data + offset;
        const double* const ax = AX.data() + offset;
        std::vector<double> x_full(pairs, 0.0);
        std::vector<double> ax_full(pairs, 0.0);
#ifdef __MPI
        LR_Util::gather_2d_to_full(px, x, x_full.data(), false, nv, no);
        LR_Util::gather_2d_to_full(px, ax, ax_full.data(), false, nv, no);
#else
        const double* const x_end = x + pairs;
        const double* const ax_end = ax + pairs;
        std::copy(x, x_end, x_full.begin());
        std::copy(ax, ax_end, ax_full.begin());
#endif
        const auto W = window_multipliers(no, nv, nactive, virtual_energies,
            x_full.data(), ax_full.data());
        std::vector<double> w_local(vv_storage, 0.0);
        std::vector<double> ew_local(w_local.size(), 0.0);
        for (int col = 0; col < pvv.ncloc(); ++col)
        {
            const int a = pvv.icglb(col);
            for (int row = 0; row < pvv.nrloc(); ++row)
            {
                const int r = pvv.irglb(row);
                const int index = col * pvv.nrloc() + row;
                w_local[index] = W[a * nv + r];
                ew_local[index] = virtual_energies[a] * w_local[index];
            }
        }
        result.density.push_back(transform_window_matrix(inputs, c, pvv, w_local.data(), LR_Util::MO_TYPE::VV, 0));
        result.direct_edm.push_back(transform_window_matrix(inputs, c, pvv, ew_local.data(), LR_Util::MO_TYPE::VV, 0));
    }
    inputs.ofs << "Window response: " << nactive << " / " << nv
        << " virtual orbitals; " << nstates << " " << inputs.spin_types[ispin]
        << " amplitudes; Dwin and direct EDM constructed; GS Hessian uses the combined density." << std::endl;
    ModuleBase::timer::end("LR", "build_window_response");
    return result;
}

template <>
WindowResponse<std::complex<double>> build_window_response(
    const GradientInputs<std::complex<double>>& inputs, const container::Tensor& X,
    const int nstates, const int nactive, const int ispin, const std::string& in_dir)
{
    if (nactive != inputs.nvirt[0])
    {
        throw std::runtime_error("Window response currently requires real Gamma orbitals");
    }
    return {};
}

template <>
WindowResponse<double> build_open_shell_window_response(const GradientInputs<double>& inputs,
    const container::Tensor& X, const int nstates, const std::vector<int>& nactive)
{
    ModuleBase::timer::start("LR", "build_open_shell_window_response");
    WindowResponse<double> result;
    if (nactive == inputs.nvirt)
    {
        ModuleBase::timer::end("LR", "build_open_shell_window_response");
        return result;
    }
    if (inputs.nspin != 2 || inputs.nk != 1 || nactive.size() != 2)
    {
        throw std::runtime_error("Open-shell window response requires two spins and real Gamma orbitals");
    }
    std::vector<psi::Psi<double>> c;
    std::vector<Parallel_Blkcyclic> pvv(2);
    std::vector<int> local_size(2);
    for (int spin = 0; spin < 2; ++spin)
    {
        c.push_back(LR_Util::get_psi_spin(inputs.psi_ks, spin, inputs.nk));
        const Parallel_Blkcyclic& px = inputs.px[spin];
        local_size[spin] = px.size();
        LR_Util::setup_2d_division(pvv[spin], px.nb(),
            inputs.nvirt[spin], inputs.nvirt[spin]
#ifdef __MPI
            , px.ictxt()
#endif
        );
    }
    const std::vector<int> spin_offset = {0, local_size[0]};
    HamiltULR<double> full_casida(inputs.xc_kernel, inputs.nspin, inputs.nbasis,
        inputs.nocc, inputs.nvirt, inputs.ucell, inputs.orb_cutoff, inputs.gd,
        inputs.psi_ks, inputs.eig_ks,
#ifdef __EXX
        inputs.exx_lri, inputs.hybrid_alpha,
#endif
        inputs.pot, inputs.kv, inputs.px, inputs.pc, inputs.pmat);
    const int batch_size = nstates * inputs.nloc;
    const int batch_storage = std::max(batch_size, 1);
    std::vector<double> AX(batch_storage, 0.0);
    const double empty_block = 0.0;
    const double* const x_data = X.data<double>() == nullptr ? &empty_block : X.data<double>();
    // Both spin blocks enter AX. Do not normalize the individual spin amplitudes.
    full_casida.hPsi(x_data, AX.data(), inputs.nloc, nstates);

    result.density.resize(nstates);
    result.direct_edm.resize(nstates);
    for (int state = 0; state < nstates; ++state)
    {
        const int state_offset = state * inputs.nloc;
        for (int spin = 0; spin < 2; ++spin)
        {
            const int no = inputs.nocc[spin];
            const int nv = inputs.nvirt[spin];
            const int pairs = no * nv;
            const int offset = state_offset + spin_offset[spin];
            const double* const x = x_data + offset;
            const double* const ax = AX.data() + offset;
            const Parallel_Blkcyclic& px = inputs.px[spin];
            std::vector<double> x_full(pairs, 0.0);
            std::vector<double> ax_full(pairs, 0.0);
#ifdef __MPI
            LR_Util::gather_2d_to_full(px, x, x_full.data(), false, nv, no);
            LR_Util::gather_2d_to_full(px, ax, ax_full.data(), false, nv, no);
#else
            const double* const x_end = x + pairs;
            const double* const ax_end = ax + pairs;
            std::copy(x, x_end, x_full.begin());
            std::copy(ax, ax_end, ax_full.begin());
#endif
            const double* const virtual_energies = inputs.eig_ks.c + spin * inputs.eig_ks.nc + no;
            const auto W = window_multipliers(no, nv, nactive[spin], virtual_energies,
                x_full.data(), ax_full.data());
            const Parallel_Blkcyclic& pv = pvv[spin];
            const int vv_storage = std::max<int>(pv.size(), 1);
            std::vector<double> w_local(vv_storage, 0.0);
            std::vector<double> ew_local(vv_storage, 0.0);
            for (int col = 0; col < pv.ncloc(); ++col)
            {
                const int a = pv.icglb(col);
                for (int row = 0; row < pv.nrloc(); ++row)
                {
                    const int r = pv.irglb(row);
                    const int index = col * pv.nrloc() + row;
                    w_local[index] = W[a * nv + r];
                    ew_local[index] = virtual_energies[a] * w_local[index];
                }
            }
            const auto density = transform_window_matrix(inputs, c[spin], pv,
                w_local.data(), LR_Util::MO_TYPE::VV, spin);
            const auto edm = transform_window_matrix(inputs, c[spin], pv,
                ew_local.data(), LR_Util::MO_TYPE::VV, spin);
            result.density[state].push_back(density[0]);
            result.direct_edm[state].push_back(edm[0]);
        }
    }
    inputs.ofs << "Open-shell window response: up " << nactive[0] << " / " << inputs.nvirt[0]
        << ", down " << nactive[1] << " / " << inputs.nvirt[1] << "; " << nstates
        << " amplitudes; spin densities and direct EDM constructed; GS Hessian uses combined densities." << std::endl;
    ModuleBase::timer::end("LR", "build_open_shell_window_response");
    return result;
}

template <>
WindowResponse<std::complex<double>> build_open_shell_window_response(
    const GradientInputs<std::complex<double>>& inputs, const container::Tensor& X,
    const int nstates, const std::vector<int>& nactive)
{
    if (nactive != inputs.nvirt)
    {
        throw std::runtime_error("Open-shell window response currently requires real Gamma orbitals");
    }
    return {};
}
}
