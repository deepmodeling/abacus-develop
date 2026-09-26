#include <vector>
#include <complex>

#include "source_pw/module_proj/onsite_proj.h"
#include "source_pw/module_proj/onsite_proj_print.h"
#include "source_cell/cell_tools.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_estate/occ_comput.h"
#include "source_io/module_parameter/parameter.h"

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::tabulate_atomic(const int ik, const char grad)
{
    ModuleBase::timer::start("OnsiteProj", "tabulate_atomic");
    // The actual tabulation of <G+k|alpha_i> (STAGE 1 + STAGE 2) is performed by
    // Onsite_Proj_tools; this member only records the k-point dimensions.
    this->ik_ = ik;
    this->npw_ = pw_basis_->npwk[ik];
    this->npwx_ = pw_basis_->npwk_max;
    ModuleBase::timer::end("OnsiteProj", "tabulate_atomic");
}

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::overlap_proj_psi(
                    const int npm,
                    const std::complex<double>* ppsi,
                    const int ld_psi)
{
    ModuleBase::timer::start("OnsiteProj", "overlap");
    // notes on refactor for DCU calculation
    // the npm here is nbands(occ) * npol, for calling cal_becp, the npol should be divided.
    int npol = this->ucell->get_npol();
    if(this->becp == nullptr || this->size_becp < npm*this->tot_nproj)
    {
        this->size_becp = npm*this->tot_nproj;
        resmem_complex_op()(this->becp, this->size_becp);
        if(this->device == base_device::GpuDevice )
        {
            resmem_complex_h_op()(this->h_becp, this->size_becp);
        }
        else
        {
            this->h_becp = this->becp;
        }
    }
    this->fs_tools->cal_becp(ik_, npm/npol, this->becp, ppsi, ld_psi > 0 ? ld_psi : this->npwx_); // in cal_becp, npm should be the one not multiplied by npol
    if(this->device == base_device::GpuDevice)
    {
        syncmem_complex_d2h_op()(h_becp, this->becp, this->size_becp);
    }
    ModuleBase::timer::end("OnsiteProj", "overlap");
}

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::cal_occupations(
        const psi::Psi<std::complex<T>, Device>* psi_in,
        const ModuleBase::matrix& wg_in,
        const int nspin_in)
{
    ModuleBase::timer::start("OnsiteProj", "cal_occupation");
    this->tabulate_atomic(0);
    std::vector<std::complex<double>> occs(this->tot_nproj * 4, 0.0);

    // loop over k-points to calculate Mi of \sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
    const int nbands = psi_in->get_nbands();
    const int npol = psi_in->get_npol();
    for(int ik = 0; ik < psi_in->get_nk(); ik++)
    {
        psi_in->fix_k(ik);
        if(ik != 0)
        {
            this->tabulate_atomic(ik);
        }
        this->overlap_proj_psi(nbands * npol, psi_in->get_pointer());
        // proj(nbands*npol , nkb) holds <alpha_{iprj}|Psi_{k,i}>.
        // nspin=2 (npol=1): the spin-up and spin-down channels are separate
        // k-points, selected by isk. nspin=1 (npol=1): no spin polarization,
        // the occupancy is split evenly so the printed magnetization is zero.
        // nspin=4 (npol=2): both spinor components are interleaved per band.
        const std::complex<double>* proj_p = this->get_h_becp();
        const double* wg_ik = &wg_in(ik, 0);
        const int isk = (nspin_in == 2 && this->isk_ != nullptr) ? this->isk_[ik] : 0;
        const int nat = static_cast<int>(this->iat_nh.size());
        elecstate::occ_from_proj(
            proj_p,
            wg_ik,
            nbands,
            npol,
            this->tot_nproj,
            nspin_in,
            isk,
            this->iat_nh.data(),
            nat,
            occs.data());
    }
    // reduce mag from all k-pools
    const int npool = GlobalV::KPAR * PARAM.inp.bndpar;
    Parallel_Reduce::reduce_double_allpool(npool, GlobalV::NPROC_IN_POOL, (double*)(&(occs[0])), occs.size()*2);
    // occ has been reduced and calculate mag
    // Print orbital charge analysis
    auto atom_labels = unitcell::get_atomLabels(this->ucell->atoms, this->ucell->ntype);
    print::print_orb_chg(this->ucell, occs, this->iat_nh, atom_labels);

    // print charge
    ModuleBase::timer::end("OnsiteProj", "cal_occupation");
}

// explicit method instantiation
template
void projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::tabulate_atomic(
    const int, const char);

template
void projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::overlap_proj_psi(
    const int, const std::complex<double>*, const int);

template
void projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::cal_occupations(
    const psi::Psi<std::complex<double>, base_device::DEVICE_CPU>*,
    const ModuleBase::matrix&,
    const int);

#if ((defined __CUDA) || (defined __ROCM))
template
void projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::tabulate_atomic(
    const int, const char);

template
void projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::overlap_proj_psi(
    const int, const std::complex<double>*, const int);

template
void projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::cal_occupations(
    const psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*,
    const ModuleBase::matrix&,
    const int);
#endif
