#include <vector>
#include <complex>

#include "source_pw/module_pwdft/onsite_proj.h"
#include "source_pw/module_pwdft/onsite_proj_print.h"
#include "source_cell/cell_tools.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::tabulate_atomic(const int ik, const char grad)
{
    ModuleBase::timer::start("OnsiteProj", "tabulate_atomic");
    // assert(grad == 'n' || grad == 'x' || grad == 'y' || grad == 'z');
    // grad = 'n' means no gradient, grad = 'x' means gradient along x, etc.

    // STAGE 1 - calculate the <G+k|p> for the given G+k vector
    // CACHE 1 - if cache the tab_, <G+k|p> can be reused for SCF and RELAX calculation
    // [in] pw_basis, ik, omega, tpiba, irow2it
    this->ik_ = ik;
    this->npw_ = pw_basis_->npwk[ik];
    this->npwx_ = pw_basis_->npwk_max;
    // std::vector<ModuleBase::Vector3<double>> q(this->npw_);
    // for(int ig = 0; ig < this->npw_; ++ig)
    // {
    //     q[ig] = pw_basis_->getgpluskcar(ik, ig); // get the G+k vector, G+k will change during CELL-RELAX
    // }
    // const int nrow = irow2it_.size();
    // std::vector<std::complex<double>> tab_(nrow*this->npw_);
    // // convention used here: 'l': <p|G+k>, 'r': <G+k|p>
    // // denote q=G+k, <r|q> = exp(iqr), the routine Fourier Transform written as F(q) = <q|f>
    // // what is calculated is <p|q> here
    // rp_.sbtft(q, tab_, 'l', this->ucell->omega, this->ucell->tpiba);

    // STAGE 2 - make_atomic: multiply e^iqtau and extend the <G+k|p> to <G+k|pi> for each atom
    // CACHE 2 - if cache the tab_atomic_, <G+k|p> can be reused for SCF calculation
    // [in] it2ia, itiaiprojm2irow, tab_, npw, sf
    // for(int irow = 0; irow < nrow; ++irow)
    // {
    //     const int it = irow2it_[irow];
    //     const int iproj = irow2iproj_[irow];
    //     const int m = irow2m_[irow];
    //     for(int ia = 0; ia < na[it]; ++ia)
    //     {
    //         // why Structure_Factor needs the FULL pw_basis???
    //         std::complex<double>* sk = this->sf_->get_sk(ik, it, ia, pw_basis_); // exp(-iqtau)
    //         // Note: idea on extending the param list of get_sk
    //         // the get_sk should have an extra param 'grad' to calculate the gradient of S(q), which
    //         // is actually very simple to be
    //         // d(S(q))/dq = -i S(q) * tau, for one direction it is just -i S(q) * tau_x (if x is the direction)
    //         const int irow_out = itiaiprojm2irow_.at(std::make_tuple(it, ia, iproj, m));
    //         for(int ig = 0; ig < this->npw_; ++ig)
    //         {
    //             std::complex<double> deriv = (grad == 'n')? 1.0: ModuleBase::NEG_IMAG_UNIT; // because sk is exp(-iqtau)
    //             deriv = (grad == 'n')? 1.0: (grad == 'x')? deriv * q[ig].x: (grad == 'y')? deriv * q[ig].y: deriv * q[ig].z;
    //             // there must be something twisted in ABACUS
    //             // because the tab_ is <p|G+k>, but the sk is exp(-iqtau). How can it get the
    //             // correct result?
    //             this->tab_atomic_[irow_out*this->npw_ + ig] = sk[ig] * tab_[irow*this->npw_ + ig] * deriv;
    //         }
    //         delete[] sk;
    //     }
    // }
    // q.clear();
    // q.shrink_to_fit();    // release memory
    // tab_.clear();
    // tab_.shrink_to_fit(); // release memory
    ModuleBase::timer::end("OnsiteProj", "tabulate_atomic");
}

template<typename T, typename Device>
void projectors::OnsiteProjector<T, Device>::overlap_proj_psi(
                    const int npm,
                    const std::complex<double>* ppsi,
                    const int ld_psi)
{
    ModuleBase::timer::start("OnsiteProj", "overlap");
    // STAGE 3 - cal_becp
    // CACHE 3 - it is no use to cache becp, it will change in each SCF iteration
    // [in] psi, tab_atomic_, npw, becp, ik
//     const char transa = 'C';
//     const char transb = 'N';
//     const int ldb = this->npwx_;
//     const int ldc = this->tot_nproj;
//     const std::complex<double> alpha = 1.0;
//     const std::complex<double> beta = 0.0;
//     if(this->becp == nullptr || this->size_becp < npm*ldc)
//     {
//         delete[] this->becp;
//         this->becp = new std::complex<double>[npm*ldc];
//         this->size_becp = npm*ldc;
//     }
//     setmem_complex_op()(ctx, this->becp, 0.0, this->size_becp);
//     gemm_op()(
//         this->ctx,
//         transa,                 // const char transa
//         transb,                 // const char transb
//         ldc,                    // const int m
//         npm,                    // const int n
//         this->npw_,             // const int k
//         &alpha,                 // const std::complex<double> alpha
//         this->tab_atomic_,      // const std::complex<double>* a
//         this->npw_,             // const int lda
//         ppsi,                   // const std::complex<double>* b
//         ldb,                    // const int ldb
//         &beta,                  // const std::complex<double> beta
//         becp,                   // std::complex<double>* c
//         ldc);                   // const int ldc
// #ifdef __MPI
//     Parallel_Reduce::reduce_pool(becp, size_becp);
// #endif

    // notes on refactor for DCU calculation
    // the npm here is nbands(occ) * npol, for calling cal_becp, the npol should be divided.
    // std::cout << "npm: " << npm << std::endl;
    // std::cout << "at " << __FILE__ << ": " << __LINE__ << " output tot_nproj: " << this->tot_nproj << std::endl;
    // std::cout << "at " << __FILE__ << ": " << __LINE__ << " output npm: " << npm << std::endl;
    // std::cout << "at " << __FILE__ << ": " << __LINE__ << " ik_: " << ik_ << std::endl;
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
        const ModuleBase::matrix& wg_in)
{
    ModuleBase::timer::start("OnsiteProj", "cal_occupation");
    this->tabulate_atomic(0);
    std::vector<std::complex<double>> occs(this->tot_nproj * 4, 0.0);

    // loop over k-points to calculate Mi of \sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
    const int nbands = psi_in->get_nbands();
    for(int ik = 0; ik < psi_in->get_nk(); ik++)
    {
        psi_in->fix_k(ik);
        if(ik != 0)
        {
            this->tabulate_atomic(ik);
        }
        // std::cout << __FILE__ << ":" << __LINE__ << " nbands = " << nbands << std::endl;
        this->overlap_proj_psi(
                        nbands * psi_in->get_npol(),
                        psi_in->get_pointer());
        const std::complex<double>* becp_p = this->get_h_becp();
        // becp(nbands*npol , nkb)
        // mag = wg * \sum_{nh}becp * becp
        int nkb = this->tot_nproj;
        //nkb = 18;
        //std::cout << "at " << __FILE__ << ": " << __LINE__ << " output nbands: " << nbands << std::endl;
        //std::cout << "at " << __FILE__ << ": " << __LINE__ << " output nkb: " << nkb << std::endl;
        for(int ib = 0;ib<nbands;ib++)
        {
            const double weight = wg_in(ik, ib);
            int begin_ih = 0;
            for(int iat = 0; iat < this->iat_nh.size(); iat++)
            {
                const int nh = this->get_nh(iat);
                for(int ih = 0; ih < nh; ih++)
                {
                    const int occ_index = (begin_ih + ih) * 4;
                    const int index = ib*2*nkb + begin_ih + ih;
                    occs[occ_index] += weight * conj(becp_p[index]) * becp_p[index];
                    occs[occ_index + 1] += weight * conj(becp_p[index]) * becp_p[index + nkb];
                    occs[occ_index + 2] += weight * conj(becp_p[index + nkb]) * becp_p[index];
                    occs[occ_index + 3] += weight * conj(becp_p[index + nkb]) * becp_p[index + nkb];
                }
                begin_ih += nh;
            }
        }
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
    const ModuleBase::matrix&);

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
    const ModuleBase::matrix&);
#endif
