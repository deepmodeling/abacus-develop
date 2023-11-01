#include "nonlocal_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include "module_base/parallel_reduce.h"
#include "module_base/tool_quit.h"
#include "module_psi/kernels/device.h"
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

namespace hamilt {

template<typename T, typename Device>
Nonlocal<OperatorPW<T, Device>>::Nonlocal(const int* isk_in,
                                               const pseudopot_cell_vnl* ppcell_in,
                                               const UnitCell* ucell_in,
                                               const ModulePW::PW_Basis_K* wfc_basis)
{
    this->classname = "Nonlocal";
    this->cal_type = pw_nonlocal;
    this->wfcpw = wfc_basis;
    this->isk = isk_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->deeq = this->ppcell->template get_deeq_data<Real>();
    this->deeq_nc = this->ppcell->template get_deeq_nc_data<Real>();
    this->vkb = this->ppcell->template get_vkb_data<Real>();
    if( this->isk == nullptr || this->ppcell == nullptr || this->ucell == nullptr)
    {
        ModuleBase::WARNING_QUIT("NonlocalPW", "Constuctor of Operator::NonlocalPW is failed, please check your code!");
    }
}

template<typename T, typename Device>
Nonlocal<OperatorPW<T, Device>>::~Nonlocal() {
    delmem_complex_op()(this->ctx, this->ps);
    delmem_complex_op()(this->ctx, this->becp);
}

template<typename T, typename Device>
void Nonlocal<OperatorPW<T, Device>>::init(const int ik_in)
{
    ModuleBase::timer::tick("Nonlocal", "getvnl");
    this->ik = ik_in;
    // Calculate nonlocal pseudopotential vkb
	if(this->ppcell->nkb > 0) //xiaohui add 2013-09-02. Attention...
	{
		this->ppcell->getvnl(this->ctx, this->ik, this->vkb);
	}

    if(this->next_op != nullptr)
    {
        this->next_op->init(ik_in);
    }

    ModuleBase::timer::tick("Nonlocal", "getvnl");
}

//--------------------------------------------------------------------------
// this function sum up each non-local pseudopotential located on each atom,
//--------------------------------------------------------------------------
template<typename T, typename Device>
void Nonlocal<OperatorPW<T, Device>>::add_nonlocal_pp(T *hpsi_in, const T *becp, const int m) const
{
    ModuleBase::timer::tick("Nonlocal", "add_nonlocal_pp");

    // number of projectors
    int nkb = this->ppcell->nkb;

    // T *ps = new T[nkb * m];
    // ModuleBase::GlobalFunc::ZEROS(ps, m * nkb);
    if (this->nkb_m < m * nkb) {
        resmem_complex_op()(this->ctx, this->ps, nkb * m, "Nonlocal<PW>::ps");
        this->nkb_m = m * nkb;
    }
    setmem_complex_op()(this->ctx, this->ps, 0, nkb * m);

    int sum = 0;
    int iat = 0;
    if (this->npol == 1)
    {
        const int current_spin = this->isk[this->ik];
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            const int nproj = this->ucell->atoms[it].ncpp.nh;
            // denghui replace 2022-10-20
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            nonlocal_op()(
                this->ctx,   // device context
                this->ucell->atoms[it].na, m, nproj, // four loop size
                sum, iat, current_spin, nkb,   // additional index params
                this->ppcell->deeq.getBound2(), this->ppcell->deeq.getBound3(), this->ppcell->deeq.getBound4(), // realArray operator()
                this->deeq, // array of data
                this->ps, this->becp); //  array of data
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            // for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            // {
            //     // each atom has nproj, means this is with structure factor;
            //     // each projector (each atom) must multiply coefficient
            //     // with all the other projectors.
            //     for (int ib = 0; ib < m; ++ib)
            //     {
            //         for (int ip2 = 0; ip2 < nproj; ip2++)
            //         {
            //             for (int ip = 0; ip < nproj; ip++)
            //             {
            //                 this->ps[(sum + ip2) * m + ib]
            //                     += this->ppcell->deeq(current_spin, iat, ip, ip2) * this->becp[ib * nkb + sum + ip];
            //             } // end ib
            //         } // end ih
            //     } // end jh
            //     sum += nproj;
            //     ++iat;
            // } // end na
        } // end nt
    }
    else
    {
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            const int nproj = this->ucell->atoms[it].ncpp.nh;
            // added by denghui at 20221109
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            nonlocal_op()(
                this->ctx,   // device context
                this->ucell->atoms[it].na, m, nproj, // four loop size
                sum, iat, nkb,   // additional index params
                this->ppcell->deeq_nc.getBound2(), this->ppcell->deeq_nc.getBound3(), this->ppcell->deeq_nc.getBound4(), // realArray operator()
                this->deeq_nc, // array of data
                this->ps, this->becp); //  array of data
            // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            // for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            // {
            //     // each atom has nproj, means this is with structure factor;
            //     // each projector (each atom) must multiply coefficient
            //     // with all the other projectors.
            //     for (int ib = 0; ib < m; ib+=2)
            //     {
            //         for (int ip2 = 0; ip2 < nproj; ip2++)
            //         {
            //             for (int ip = 0; ip < nproj; ip++)
            //             {
            //                 psind = (sum + ip2) * m + ib;
            //                 becpind = ib * nkb + sum + ip;
            //                 becp1 = becp[becpind];
            //                 becp2 = becp[becpind + nkb];
            //                 ps[psind] += this->ppcell->deeq_nc(0, iat, ip2, ip) * becp1
            //                              + this->ppcell->deeq_nc(1, iat, ip2, ip) * becp2;
            //                 ps[psind + 1] += this->ppcell->deeq_nc(2, iat, ip2, ip) * becp1
            //                                  + this->ppcell->deeq_nc(3, iat, ip2, ip) * becp2;
            //             } // end ib
            //         } // end ih
            //     } // end jh
            //     sum += nproj;
            //     ++iat;
            // } // end na
        } // end nt
    }

    // use simple method.
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui optimize 2021-3-31
    char transa = 'N';
    char transb = 'T';
    if (m == 1)
    {
        int inc = 1;
        // denghui replace 2022-10-20
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        gemv_op()(
            this->ctx,
            transa,
            this->npw,
            this->ppcell->nkb,
            &this->one,
            this->vkb,
            this->ppcell->vkb.nc,
            this->ps,
            inc,
            &this->one,
            hpsi_in,
            inc);
    }
    else
    {
        int npm = m;
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // denghui replace 2022-10-20
        gemm_op()(
            this->ctx,
            transa,
            transb,
            this->npw,
            npm,
            this->ppcell->nkb,
            &this->one,
            this->vkb,
            this->ppcell->vkb.nc,
            this->ps,
            npm,
            &this->one,
            hpsi_in,
            this->max_npw
        );
    }
    ModuleBase::timer::tick("Nonlocal", "add_nonlocal_pp");
}

template<typename T, typename Device>
void Nonlocal<OperatorPW<T, Device>>::act(
    const int nbands,
    const int nbasis,
    const int npol,
    const T* tmpsi_in,
    T* tmhpsi,
    const int ngk_ik)const
{
    ModuleBase::timer::tick("Operator", "NonlocalPW");
    if(!GlobalV::use_paw)
    {
        this->npw = ngk_ik;
        this->max_npw = nbasis / npol;
        this->npol = npol;

        if (this->ppcell->nkb > 0)
        {
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            // qianrui optimize 2021-3-31
            int nkb = this->ppcell->nkb;
            if (this->nkb_m < nbands * nkb) {
                resmem_complex_op()(this->ctx, this->becp, nbands * nkb, "Nonlocal<PW>::becp");
            }
            // ModuleBase::ComplexMatrix becp(nbands, nkb, false);
            char transa = 'C';
            char transb = 'N';
            if (nbands == 1)
            {
                int inc = 1;
                // denghui replace 2022-10-20
                // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                gemv_op()(
                    this->ctx,
                    transa,
                    this->npw,
                    nkb,
                    &this->one,
                    this->vkb,
                    this->ppcell->vkb.nc,
                    tmpsi_in,
                    inc,
                    &this->zero,
                    this->becp,
                    inc);
            }
            else
            {
                int npm = nbands;
                //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                // denghui replace 2022-10-20
                gemm_op()(
                    this->ctx,
                    transa,
                    transb,
                    nkb,
                    npm,
                    this->npw,
                    &this->one,
                    this->vkb,
                    this->ppcell->vkb.nc,
                    tmpsi_in,
                    max_npw,
                    &this->zero,
                    this->becp,
                    nkb
                );
            }

            Parallel_Reduce::reduce_pool(becp, nkb * nbands);

            this->add_nonlocal_pp(tmhpsi, becp, nbands);
        }
    }
    else
    {
#ifdef USE_PAW
        this->npw = ngk_ik;
        this->max_npw = nbasis / npol;
        this->npol = npol;
        std::complex<double> *vnlpsi;
        vnlpsi = new std::complex<double> [npw];
        for(int ibands = 0; ibands < nbands; ibands++)
        {
            GlobalC::paw_cell.paw_nl_psi(0,reinterpret_cast<const std::complex<double>*> (&tmpsi_in[ibands*max_npw]),vnlpsi);
            for(int i = 0; i < npw; i++)
            {
                tmhpsi[ibands*max_npw+i] += vnlpsi[i];
            }
        }
#endif
    }
    ModuleBase::timer::tick("Operator", "NonlocalPW");
}

template<typename T, typename Device>
template<typename T_in, typename Device_in>
hamilt::Nonlocal<OperatorPW<T, Device>>::Nonlocal(const Nonlocal<OperatorPW<T_in, Device_in>> *nonlocal)
{
    this->classname = "Nonlocal";
    this->cal_type = pw_nonlocal;
    this->ik = nonlocal->get_ik();
    this->isk = nonlocal->get_isk();
    this->ppcell = nonlocal->get_ppcell();
    this->ucell = nonlocal->get_ucell();
    this->deeq = this->ppcell->d_deeq;
    this->deeq_nc = this->ppcell->template get_deeq_nc_data<Real>();
    this->vkb = this->ppcell->template get_vkb_data<Real>();
    if( this->isk == nullptr || this->ppcell == nullptr || this->ucell == nullptr)
    {
        ModuleBase::WARNING_QUIT("NonlocalPW", "Constuctor of Operator::NonlocalPW is failed, please check your code!");
    }
}

template class Nonlocal<OperatorPW<std::complex<float>, psi::DEVICE_CPU>>;
template class Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_CPU>>;
// template Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_CPU>>::Nonlocal(const Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_CPU>> *nonlocal);
#if ((defined __CUDA) || (defined __ROCM))
template class Nonlocal<OperatorPW<std::complex<float>, psi::DEVICE_GPU>>;
template class Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_GPU>>;
// template Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_CPU>>::Nonlocal(const Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_GPU>> *nonlocal);
// template Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_GPU>>::Nonlocal(const Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_CPU>> *nonlocal);
// template Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_GPU>>::Nonlocal(const Nonlocal<OperatorPW<std::complex<double>, psi::DEVICE_GPU>> *nonlocal);
#endif
} // namespace hamilt