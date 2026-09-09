#include "source_estate/uspp_density.h"

#include "source_base/constants.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/libm/libm.h"
#include "source_base/macros.h"
#include "source_base/math_ylmreal.h"
#include "source_base/module_container/ATen/core/tensor.h"
#include "source_base/parallel_reduce.h"

#include <cassert>

namespace elecstate
{
int uspp_becsum_size(const UnitCell& ucell, const pseudopot_cell_vnl& ppcell, const int nspin)
{
    // Store only i <= j projector pairs, with the same maximum-sized block for every atom.
    return nspin * ucell.nat * ppcell.nhm * (ppcell.nhm + 1) / 2;
}

// The implementation owns device tensors while exposing only host coefficients to callers.
template <typename T, typename Device>
class UsppProjector<T, Device>::Workspace
{
  public:
    using ContainerDevice = typename ct::PsiToContainer<Device>::type;
    const ct::DeviceType device = ct::DeviceTypeToEnum<ContainerDevice>::value;
    ct::Tensor overlaps;
    ct::Tensor left;
    ct::Tensor right;
    ct::Tensor products;
    std::vector<T> overlaps_host;
    std::vector<T> left_host;
    std::vector<T> right_host;
    std::vector<T> products_host;

    Workspace(const int nbands, const int nkb, const int nhm)
        : overlaps(ct::DataTypeToEnum<T>::value, device, ct::TensorShape({nbands * nkb})),
          left(ct::DataTypeToEnum<T>::value, device, ct::TensorShape({nbands * nhm})),
          right(ct::DataTypeToEnum<T>::value, device, ct::TensorShape({nbands * nhm})),
          products(ct::DataTypeToEnum<T>::value, device, ct::TensorShape({nhm * nhm})), overlaps_host(nbands * nkb),
          left_host(nbands * nhm), right_host(nbands * nhm), products_host(nhm * nhm)
    {
    }
};

template <typename T, typename Device>
UsppProjector<T, Device>::UsppProjector(const UnitCell& ucell, const pseudopot_cell_vnl& ppcell, const int nbands)
    : ucell_(ucell), ppcell_(ppcell), nbands_(nbands), work_(new Workspace(nbands, ppcell.nkb, ppcell.nhm))
{
}

// Defining this destructor here makes Workspace complete when unique_ptr releases it.
template <typename T, typename Device>
UsppProjector<T, Device>::~UsppProjector() = default;

template <typename T, typename Device>
void UsppProjector<T, Device>::accumulate(const int ik,
                                          const T* coefficients,
                                          const int stride,
                                          const int npw,
                                          const int spin,
                                          const std::vector<double>& weights,
                                          std::vector<double>* becsum)
{
    using Real = typename GetTypeReal<T>::type;
    // Memory-copy template arguments specify the destination device first, then the source.
    using CopyToHost = base_device::memory::synchronize_memory_op<T, base_device::DEVICE_CPU, Device>;
    using CopyToDevice = base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>;
    assert(static_cast<int>(weights.size()) == nbands_);
    if (ppcell_.nkb == 0 || nbands_ == 0)
    {
        return;
    }
    Workspace& work = *work_;
    const T one(1, 0);
    const T zero(0, 0);
    const int nkb = ppcell_.nkb;
    const int nh_tot = ppcell_.nhm * (ppcell_.nhm + 1) / 2;
    assert(spin >= 0 && static_cast<int>(becsum->size()) >= (spin + 1) * ucell_.nat * nh_tot);
    // Build beta_i(G+k) for this k point before projecting the pseudo-wavefunctions.
    T* vkb = ppcell_.template get_vkb_data<Real>();
    ppcell_.getvnl(static_cast<Device*>(nullptr), ucell_, ik, vkb);
    // Each rank forms its PW contribution to <beta|psi>; the pool sum completes the overlap.
    // BLAS 'C' conjugate-transposes the projector matrix; stride includes any padding between states.
    T* overlaps = work.overlaps.template data<T>();
    ModuleBase::gemm_op<T, Device>()('C', 'N', nkb, nbands_, npw, &one, vkb, ppcell_.vkbnc, coefficients, stride, &zero, overlaps, nkb);
    CopyToHost()(work.overlaps_host.data(), overlaps, nbands_ * nkb);
    Parallel_Reduce::reduce_pool(work.overlaps_host.data(), nbands_ * nkb);

    for (int it = 0; it < ucell_.ntype; ++it)
    {
        const Atom& atom = ucell_.atoms[it];
        if (!atom.ncpp.tvanp)
        {
            continue;
        }
        const int nh = atom.ncpp.nh;
        for (int ia = 0; ia < atom.na; ++ia)
        {
            const int iat = ucell_.itia2iat(it, ia);
            // Gather this atom's overlaps into band-by-projector matrices for the band sum below.
            // Callers supply the state weights: occupations and k weights for SCF, selected-state weights for output.
            for (int ih = 0; ih < nh; ++ih)
            {
                const int ikb = ppcell_.indv_ijkb0[iat] + ih;
                for (int ib = 0; ib < nbands_; ++ib)
                {
                    work.left_host[ih * nbands_ + ib] = work.overlaps_host[ib * nkb + ikb];
                    work.right_host[ih * nbands_ + ib] = work.left_host[ih * nbands_ + ib] * static_cast<Real>(weights[ib]);
                }
            }
            // The small matrix is sum_n weight_n conj(<beta_i|psi_n>) <beta_j|psi_n>.
            T* left = work.left.template data<T>();
            T* right = work.right.template data<T>();
            T* products = work.products.template data<T>();
            CopyToDevice()(left, work.left_host.data(), nbands_ * nh);
            CopyToDevice()(right, work.right_host.data(), nbands_ * nh);
            ModuleBase::gemm_op<T, Device>()('C', 'N', nh, nh, nbands_, &one, left, nbands_, right, nbands_, &zero, products, nh);
            CopyToHost()(work.products_host.data(), products, nh * nh);
            // Accumulate into [spin][atom][packed pair], retaining contributions from earlier calls.
            const int offset = (spin * ucell_.nat + iat) * nh_tot;
            int pair = 0;
            for (int ih = 0; ih < nh; ++ih)
            {
                for (int jh = ih; jh < nh; ++jh)
                {
                    // Q_ij(r) = Q_ji(r) is real, so conjugate off-diagonal terms give twice their real part.
                    // BLAS stores products column-major; this transposed access has the same real part.
                    const double multiplicity = ih == jh ? 1.0 : 2.0;
                    (*becsum)[offset + pair] += multiplicity * static_cast<double>(std::real(work.products_host[ih * nh + jh]));
                    ++pair;
                }
            }
        }
    }
}

void add_uspp_density(const UnitCell& ucell,
                      const pseudopot_cell_vnl& ppcell,
                      const ModulePW::PW_Basis& basis,
                      const int nspin,
                      const std::vector<double>& becsum,
                      std::complex<double>** rhog)
{

    const int npw = basis.npw;
    assert(static_cast<int>(becsum.size()) == uspp_becsum_size(ucell, ppcell, nspin));
    if (npw == 0)
    {
        return;
    }
    const int lmaxq = ppcell.lmaxq;
    const int nh_tot = ppcell.nhm * (ppcell.nhm + 1) / 2;
    const std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;

    // Add rho_aug(G) = sum_{I,i<=j} Q_ij(G) exp(-i G.R_I) B_ij^I to the existing density.
    // Q(G) already includes 1/omega; becsum holds B_ij^I, including the off-diagonal factor of two.
    std::vector<double> qmod_host(npw);
    std::vector<std::complex<double>> qgm_host(npw);
    for (int ig = 0; ig < npw; ig++)
    {
        // Convert the dimensionless reciprocal-vector length to inverse bohr for the radial Q tables.
        qmod_host[ig] = static_cast<double>(basis.gcar[ig].norm() * ucell.tpiba);
    }

    // Real spherical harmonics supply the angular dependence of Q_ij(G); radial tables supply its |G| dependence.
    ModuleBase::matrix ylmk0(lmaxq * lmaxq, npw);
    ModuleBase::YlmReal::Ylm_Real(lmaxq * lmaxq, npw, basis.gcar, ylmk0);

    for (int it = 0; it < ucell.ntype; it++)
    {
        const Atom* atom = &ucell.atoms[it];
        if (atom->ncpp.tvanp)
        {
            const int nij = atom->ncpp.nh * (atom->ncpp.nh + 1) / 2;

            // Atoms of one species share Q_ij(G); only their positions and projector coefficients differ.
            std::vector<std::complex<double>> skk_host(atom->na * npw);
            std::vector<std::complex<double>> tbecsum_host(nspin * atom->na * nij);
            for (int ia = 0; ia < atom->na; ia++)
            {
                const int iat = ucell.itia2iat(it, ia);
                for (int is = 0; is < nspin; is++)
                {
                    for (int ij = 0; ij < nij; ij++)
                    {
                        tbecsum_host[is * atom->na * nij + ia * nij + ij]
                            = static_cast<std::complex<double>>(becsum[is * ucell.nat * nh_tot + iat * nh_tot + ij]);
                    }
                }
                for (int ig = 0; ig < npw; ig++)
                {
                    // Translate the atom-centered augmentation to R_I via exp(-i G.R_I).
                    // gcar and tau use lattice-scaled units, so the phase needs the explicit factor 2*pi.
                    double arg = basis.gcar[ig] * atom->tau[ia];
                    skk_host[ia * npw + ig] = ModuleBase::libm::exp(ci_tpi * arg);
                }
            }

            for (int is = 0; is < nspin; is++)
            {
                // Sum atoms of this species: aux2(G,ij) = sum_I exp(-i G.R_I) B_ij^I.
                // In column-major BLAS, skk is [G,atom] and tbecsum is [pair,atom]; 'T' aligns the atom indices.
                std::vector<std::complex<double>> aux2_host(nij * npw);
                const std::complex<double> one_d(1, 0);
                const std::complex<double> zero_d(0, 0);
                char transa = 'N';
                char transb = 'T';
                zgemm_(&transa,
                       &transb,
                       &npw,
                       &nij,
                       &atom->na,
                       &one_d,
                       skk_host.data(),
                       &npw,
                       &tbecsum_host[is * atom->na * nij],
                       &nij,
                       &zero_d,
                       aux2_host.data(),
                       &npw);

                int ijh = 0;
                for (int ih = 0; ih < atom->ncpp.nh; ih++)
                {
                    for (int jh = ih; jh < atom->ncpp.nh; jh++)
                    {
                        // Reconstruct the atom-centered Q_ij(G) from radial tables and spherical harmonics.
                        // Use the same packed-pair order as becsum; its off-diagonal factor is already included.
                        ppcell.radial_fft_q<double,
                                            base_device::DEVICE_CPU>(nullptr, npw, ih, jh, it, qmod_host.data(), ylmk0.c, qgm_host.data());
                        for (int ig = 0; ig < npw; ig++)
                        {
                            rhog[is][ig] += qgm_host[ig] * aux2_host[ijh * npw + ig];
                        }
                        ijh++;
                    }
                }
            }
        }
    }
}

// Emit all precision/device combinations used by SCF and double-precision output.
template class UsppProjector<std::complex<float>, base_device::DEVICE_CPU>;
template class UsppProjector<std::complex<double>, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class UsppProjector<std::complex<float>, base_device::DEVICE_GPU>;
template class UsppProjector<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace elecstate
