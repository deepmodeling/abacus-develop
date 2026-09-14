#include "onsite_proj_tools.h"

#include "source_base/math_polyint.h"
#include "source_base/math_ylmreal.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_io/module_parameter/parameter.h"
#include "nonlocal_maths.hpp"

// cal_becp
// starts from vkb (nkb, ng) table
// it should be merely the multiplication of matrix (vkb, ng) * (ng, nbands) -> (vkb, nbands)
// should be irrelevant with what the matrix is.
// the vkb index generation should be maintained elsewhere.
// vkb already has atomic position information, calculated from the vq and sk
// . the multiplication with sk should be within specific operator
// because the atom selection task is operator-specific.
template <typename FPTYPE, typename Device>
void hamilt::Onsite_Proj_tools<FPTYPE, Device>::cal_becp(int ik,
                                                 int npm,
                                                 std::complex<FPTYPE>* becp_in,
                                                 const std::complex<FPTYPE>* ppsi_in,
                                                 int npwx)
{
    ModuleBase::TITLE("Onsite_Proj_tools", "cal_becp");
    ModuleBase::timer::start("Onsite_Proj_tools", "cal_becp");

    const int npol = this->ucell_->get_npol();
    const std::complex<FPTYPE>* ppsi = ppsi_in == nullptr ? &(this->psi_[0](ik, 0, 0)) : ppsi_in;
    const int npw = this->wfc_basis_->npwk[ik];
    if (becp_in == nullptr && this->becp == nullptr)
    {
        resmem_complex_op()(becp, this->nbands * npol * this->nkb);
    }
    std::complex<FPTYPE>* becp_tmp = becp_in == nullptr ? this->becp : becp_in;
    const int size_becp_act = npm * npol * this->nkb;
    if (ik != this->current_ik) // different ik, need to recalculate vkb
    {
        const int size_becp = this->nbands * npol * this->nkb;
        if (this->becp == nullptr)
        {
            resmem_complex_op()(becp, size_becp);
        }

        // prepare math tools
        Nonlocal_maths<FPTYPE, Device> maths(*(this->nhtol), this->lprojmax, this->ucell_);

        std::complex<FPTYPE>* vkb_ptr = this->ppcell_vkb;

        // calculate G+K
        this->g_plus_k = maths.cal_gk(ik, this->wfc_basis_);
        FPTYPE *gk = g_plus_k.data(), *vq_tb = this->tabtpr->ptr;
        // vq_tb has dimension (ntype, nproj, GlobalV::NQX)

        // calculate sk
        resmem_complex_op()(hd_sk, this->ucell_->nat * npw);
        this->sf_->get_sk(ctx, ik, this->wfc_basis_, hd_sk);
        std::complex<FPTYPE>* d_sk = this->hd_sk;
        // prepare ylm，size: (lmax+1)^2 * this->max_npw
        const int lmax_ = this->lprojmax;
        maths.cal_ylm(lmax_, npw, g_plus_k.data(), hd_ylm);

        // DEBUG: ONCE YOU CHECK ylm VALUES, YOU UNCOMMENT THE FOLLOW
        // std::vector<ModuleBase::Vector3<double>> qs(npw);
        // for (int ig = 0; ig < npw; ig++)
        // {
        //     qs[ig] = this->wfc_basis_->getgpluskcar(ik, ig);
        // }
        // const int total_lm = (lmax_ + 1) * (lmax_ + 1);
        // ModuleBase::matrix ylmref(total_lm, npw);
        // ModuleBase::YlmReal::Ylm_Real(total_lm, npw, qs.data(), ylmref);
        // std::cout << "Compare the Ylm values of two methods:" << std::endl;
        // int lm = 0;
        // for(int l_ = 0; l_ < lmax_ + 1; l_++)
        // {
        //     for(int m_ = -l_; m_ <= l_; m_++)
        //     {
        //         std::cout << "l = " << l_ << " m = " << m_ << std::endl;
        //         lm = l_ * l_ + l_ + m_;
        //         for(int ig = 0; ig < npw; ig++)
        //         {
        //             std::cout << "[" << ylmref(lm, ig) << " " << hd_ylm[lm * npw + ig] << "]" << " ";
        //         }
        //         std::cout << std::endl;
        //     }
        //     std::cout << std::endl;
        // }
        // ModuleBase::WARNING_QUIT("Onsite_Proj_tools", "cal_becp");

        if (this->device == base_device::GpuDevice)
        {
            syncmem_var_h2d_op()(d_g_plus_k, g_plus_k.data(), g_plus_k.size());
            syncmem_var_h2d_op()(d_vq_tab, this->tabtpr->ptr, this->tabtpr->getSize());
            gk = d_g_plus_k;
            vq_tb = d_vq_tab;
        }

        // int vkb_size = 0;
        for (int it = 0; it < this->ucell_->ntype; it++) // loop all elements
        {
            // interpolate (it, 0..nproj[it], 0..npw) to get hd_vq
            cal_vq_op()(this->ctx,
                        vq_tb, // its data is correct, dimension (ntype, nprojmax, GlobalV::NQX)
                        it,    // but maybe it is (ntype, nprojmax*npol, GlobalV::NQX)
                        gk,
                        npw,
                        this->tabtpr->getBound2(),
                        this->tabtpr->getBound3(),
                        PARAM.globalv.dq,
                        nproj[it],
                        hd_vq); // hd_vq has dimension (nprojmax, npwx), this size will be the largest needed.

            // DEBUG: ONCE YOU CHECK vq VALUES, YOU UNCOMMENT THE FOLLOWING LINE
            // for(int ip = 0; ip < nproj[it]; ip++)
            // {
            //     std::cout << "projector #" << ip << " of atom type " << it << std::endl;
            //     for(int iq = 0; iq < npw; iq++)
            //     {
            //         std::cout << hd_vq[ip * npw + iq] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << std::endl;

            // prepare（-i）^l, size: nh
            std::vector<std::complex<double>> pref = maths.cal_pref(it, h_atom_nh[it]);
            const int nh = pref.size();
            this->dvkb_indexes.resize(nh * 4);
            // print the value of nhtol
            // nhtol->print(std::cout); // as checked, nhtol works as expected
            maths.cal_dvkb_index(nproj[it], this->nhtol->c, this->nhtol->nc, npw, it, 0, 0, this->dvkb_indexes.data());

            if (this->device == base_device::GpuDevice)
            {
                syncmem_int_h2d_op()(d_dvkb_indexes, dvkb_indexes.data(), nh * 4);
                syncmem_complex_h2d_op()(d_pref_in, pref.data(), nh);
            }

            for (int ia = 0; ia < h_atom_na[it]; ia++)
            {
                if (this->device == base_device::CpuDevice)
                {
                    d_pref_in = pref.data();
                    d_dvkb_indexes = dvkb_indexes.data();
                }
                cal_vkb_op()(this->ctx, nh, npw, d_dvkb_indexes, hd_vq, hd_ylm, d_sk, d_pref_in, vkb_ptr);
                vkb_ptr += nh * npw; // vkb_ptr has dimension (nhtot, npwx), this size will be the largest needed.
                d_sk += npw;
                // vkb_size += nh * npw;
            }
        }
        this->current_ik = ik;
    }
    // DEBUG: ONCE YOU CHECK vkb VALUES, YOU UNCOMMENT THE FOLLOWING LINE
    // for (int i = 0; i < vkb_size; i++)
    // {
    //     if (i % npw == 0)
    //     {
    //         std::cout << "The #" << i / npw << " projector" << std::endl;
    //     }
    //     std::cout << this->ppcell_vkb[i] << " ";
    // }
    // std::cout << std::endl;
    // ModuleBase::WARNING_QUIT("Onsite_Proj_tools", "cal_becp");

    // PLAN: seperate the lower and upper into two parts, individually called.
    const char transa = 'C';
    const char transb = 'N';
    int npm_npol = npm * npol;
    gemm_op()(transa,
              transb,
              this->nkb,
              npm_npol, // nbands(occ)*npol
              npw,
              &ModuleBase::ONE,
              this->ppcell_vkb,
              npw,
              ppsi,
              npwx > 0 ? npwx : this->max_npw,
              &ModuleBase::ZERO,
              becp_tmp,
              this->nkb);

    if (this->device == base_device::GpuDevice)
    {
        std::complex<FPTYPE>* h_becp = nullptr;
        resmem_complex_h_op()(h_becp, size_becp_act);
        syncmem_complex_d2h_op()(h_becp, becp_tmp, size_becp_act);
        Parallel_Reduce::reduce_pool(h_becp, size_becp_act);
        syncmem_complex_h2d_op()(becp_tmp, h_becp, size_becp_act);
        delmem_complex_h_op()(h_becp);
    }
    else
    {
        Parallel_Reduce::reduce_pool(becp_tmp, size_becp_act);
    }
    // DEBUG: ONCE YOU CHECK becp VALUES, YOU UNCOMMENT THE FOLLOWING LINE
    // std::cout << "ik: " << ik << std::endl;
    // for (int i = 0; i < npm_npol*this->nkb; i++)
    // {
    //     std::cout << "becp[" << i << "]: " << becp[i] << std::endl;
    // }
    ModuleBase::timer::end("Onsite_Proj_tools", "cal_becp");
}

// explicit method instantiation
template
void hamilt::Onsite_Proj_tools<double, base_device::DEVICE_CPU>::cal_becp(
    int, int, std::complex<double>*, const std::complex<double>*, int);

#if ((defined __CUDA) || (defined __ROCM))
template
void hamilt::Onsite_Proj_tools<double, base_device::DEVICE_GPU>::cal_becp(
    int, int, std::complex<double>*, const std::complex<double>*, int);
#endif
