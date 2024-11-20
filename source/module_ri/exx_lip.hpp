//==========================================================
// AUTHOR : Peize Lin

// DATE : 2015-03-10
//==========================================================

#ifndef EXX_LIP_HPP
#define EXX_LIP_HPP

#include "exx_lip.h"
#include "module_base/vector3.h"
#include "module_base/global_function.h"
#include "module_base/vector3.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_cell/klist.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_base/lapack_connector.h"
#include "module_base/parallel_global.h"
#include "module_parameter/parameter.h"
#include "module_elecstate/elecstate.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_hamilt_pw/hamilt_pwdft/wfinit.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"

#include <limits>

template <typename T, typename Device>
void Exx_Lip<T, Device>::cal_exx()
{
    ModuleBase::TITLE("Exx_Lip", "cal_exx");
    ModuleBase::timer::tick("Exx_Lip", "cal_exx");

    this->wf_wg_cal();
    this->psi_cal();
    for (int ik = 0; ik < this->k_pack->kv_ptr->get_nks(); ++ik)
    {
        this->phi_cal(this->k_pack, ik);

        this->judge_singularity(ik);
        for (int iw_l = 0; iw_l < PARAM.globalv.nlocal; ++iw_l) {
            for (int iw_r = 0; iw_r < PARAM.globalv.nlocal; ++iw_r) {
                this->sum1[iw_l * PARAM.globalv.nlocal + iw_r] = T(0.0);
        } }
        if (Conv_Coulomb_Pot_K::Ccp_Type::Ccp == info.ccp_type || Conv_Coulomb_Pot_K::Ccp_Type::Hf == info.ccp_type)
        {
            this->sum2_factor = 0.0;
            if (gzero_rank_in_pool == GlobalV::RANK_IN_POOL) {
                for (int iw_l = 0; iw_l < PARAM.globalv.nlocal; ++iw_l) {
                    for (int iw_r = 0; iw_r < PARAM.globalv.nlocal; ++iw_r) {
                        this->sum3[iw_l][iw_r] = T(0.0);
            } } }
        }

        for (int iq_tmp = this->iq_vecik; iq_tmp < this->iq_vecik + this->q_pack->kv_ptr->get_nks() / PARAM.inp.nspin; ++iq_tmp)					// !!! k_point
            // parallel incompleted.need to loop iq in other pool
        {
            const int iq =
                (ik < (this->k_pack->kv_ptr->get_nks() / PARAM.inp.nspin))
                ? (iq_tmp % (this->q_pack->kv_ptr->get_nks() / PARAM.inp.nspin))
                : (iq_tmp % (this->q_pack->kv_ptr->get_nks() / PARAM.inp.nspin) + (this->q_pack->kv_ptr->get_nks() / PARAM.inp.nspin));
            this->qkg2_exp(ik, iq);
            for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
            {
                this->b_cal(ik, iq, ib);
                if (Conv_Coulomb_Pot_K::Ccp_Type::Ccp == info.ccp_type || Conv_Coulomb_Pot_K::Ccp_Type::Hf == info.ccp_type) {
                    if (iq == this->iq_vecik) {
                        this->sum3_cal(iq, ib);
                } }
                this->b_sum(iq, ib);
            }
        }
        this->sum_all(ik);
    }
    this->exx_energy_cal();

    ModuleBase::timer::tick("Exx_Lip", "cal_exx");
}

template <typename T, typename Device>
Exx_Lip<T, Device>::Exx_Lip(
    const Exx_Info::Exx_Info_Lip& info_in,
    const ModuleSymmetry::Symmetry& symm,
    K_Vectors* kv_ptr_in,
    psi::WFInit<T, Device>* wf_ptr_in,
    psi::Psi<T, Device>* kspw_psi_ptr_in,
    //    wavefunc* wf_ptr_in,
    const ModulePW::PW_Basis_K* wfc_basis_in,
    const ModulePW::PW_Basis* rho_basis_in,
    const Structure_Factor& sf,
    const UnitCell* ucell_ptr_in,
    const elecstate::ElecState* pelec_in) : info(info_in)
{
    ModuleBase::TITLE("Exx_Lip", "init");
    ModuleBase::timer::tick("Exx_Lip", "init");

    this->k_pack = new k_package;
    this->k_pack->kv_ptr = kv_ptr_in;
    this->k_pack->wf_ptr = wf_ptr_in;
    this->k_pack->pelec = pelec_in;
    this->k_pack->kspw_psi_ptr = kspw_psi_ptr_in;
    this->wfc_basis = wfc_basis_in;
    this->rho_basis = rho_basis_in;
    this->ucell_ptr = ucell_ptr_in;

    int gzero_judge = -1;
    if (this->rho_basis->gg_uniq[0] < 1e-8)
        { gzero_judge = GlobalV::RANK_IN_POOL; }
  #ifdef __MPI
    MPI_Allreduce(&gzero_judge, &gzero_rank_in_pool, 1, MPI_INT, MPI_MAX, POOL_WORLD);
  #endif
    this->k_pack->wf_wg.create(this->k_pack->kv_ptr->get_nks(),PARAM.inp.nbands);

    this->k_pack->hvec_array = new psi::Psi<T, Device>(this->k_pack->kv_ptr->get_nks(), PARAM.inp.nbands, PARAM.globalv.nlocal);
    // this->k_pack->hvec_array = new ModuleBase::ComplexMatrix[this->k_pack->kv_ptr->get_nks()];
    // for( int ik=0; ik<this->k_pack->kv_ptr->get_nks(); ++ik)
    // {
    // 	this->k_pack->hvec_array[ik].create(PARAM.globalv.nlocal,PARAM.inp.nbands);
    // }

    // if (PARAM.inp.init_chg=="atomic")
    {
        this->q_pack = this->k_pack;
    }
    // else if(PARAM.inp.init_chg=="file")
    // {
    //     read_q_pack(symm, this->wfc_basis, sf);
    // }

    this->phi.resize(PARAM.globalv.nlocal);
    for (int iw = 0; iw < PARAM.globalv.nlocal; ++iw)
        { this->phi[iw].resize(this->rho_basis->nrxx); }

    this->psi.resize(this->q_pack->kv_ptr->get_nks());
    for (int iq = 0; iq < this->q_pack->kv_ptr->get_nks(); ++iq)
    {
        this->psi[iq].resize(PARAM.inp.nbands);
        for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
            { this->psi[iq][ib].resize(this->rho_basis->nrxx); }
    }

    this->recip_qkg2.resize(this->rho_basis->npw);

    this->b.resize(PARAM.globalv.nlocal * this->rho_basis->npw);

    this->sum1.resize(PARAM.globalv.nlocal * PARAM.globalv.nlocal);

    if (Conv_Coulomb_Pot_K::Ccp_Type::Ccp == info.ccp_type || Conv_Coulomb_Pot_K::Ccp_Type::Hf == info.ccp_type)
    {
        if (gzero_rank_in_pool == GlobalV::RANK_IN_POOL)
        {
            this->b0.resize(PARAM.globalv.nlocal);
            this->sum3.resize(PARAM.globalv.nlocal);
            for (int iw_l = 0; iw_l < PARAM.globalv.nlocal; ++iw_l)
                { this->sum3[iw_l].resize(PARAM.globalv.nlocal); }
        }
    }

    this->exx_matrix.resize(this->k_pack->kv_ptr->get_nks());
    for (int ik = 0; ik < this->k_pack->kv_ptr->get_nks(); ++ik)
    {
        this->exx_matrix[ik].resize(PARAM.globalv.nlocal);
        for (int iw_l = 0; iw_l < PARAM.globalv.nlocal; ++iw_l)
            { this->exx_matrix[ik][iw_l].resize(PARAM.globalv.nlocal); }
    }

    ModuleBase::timer::tick("Exx_Lip", "init");
}

template <typename T, typename Device>
Exx_Lip<T, Device>::~Exx_Lip()
{
    if (this->k_pack)
        { delete this->k_pack->hvec_array;  this->k_pack->hvec_array = nullptr; }
    delete this->k_pack;	this->k_pack = nullptr;

    if (PARAM.inp.init_chg == "atomic")
    {
        this->q_pack = nullptr;
    }
    else if (PARAM.inp.init_chg == "file")
    {
        delete this->q_pack->kv_ptr;	this->q_pack->kv_ptr = nullptr;
        delete this->q_pack->wf_ptr;	this->q_pack->wf_ptr = nullptr;
        // delete[] this->q_pack->hvec_array;	this->q_pack->hvec_array=nullptr;
        delete this->q_pack;	this->q_pack = nullptr;
    }
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::wf_wg_cal()
{
    ModuleBase::TITLE("Exx_Lip", "wf_wg_cal");
    ModuleBase::timer::tick("Exx_Lip", "wf_wg_cal");
    if (PARAM.inp.nspin == 1) {
        for (int ik = 0; ik < this->k_pack->kv_ptr->get_nks(); ++ik) {
            for (int ib = 0; ib < PARAM.inp.nbands; ++ib) {
                this->k_pack->wf_wg(ik, ib) = this->k_pack->pelec->wg(ik, ib) / 2;
    } } }
    else if (PARAM.inp.nspin == 2) {
        for (int ik = 0; ik < this->k_pack->kv_ptr->get_nks(); ++ik) {
            for (int ib = 0; ib < PARAM.inp.nbands; ++ib) {
                this->k_pack->wf_wg(ik, ib) = this->k_pack->pelec->wg(ik, ib);
    } } }
    ModuleBase::timer::tick("Exx_Lip", "wf_wg_cal");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::phi_cal(k_package* kq_pack, const int ikq)
{
    ModuleBase::timer::tick("Exx_Lip", "phi_cal");
    std::vector<T> porter (this->wfc_basis->nrxx);
    for (int iw = 0; iw < PARAM.globalv.nlocal; ++iw)
    {
        // this->wfc_basis->recip2real(&kq_pack->wf_ptr->wanf2[ikq](iw,0), porter.data(), ikq);
        this->wfc_basis->recip2real(&(kq_pack->wf_ptr->get_psig().lock()->operator()(ikq, iw, 0)), porter.data(), ikq);
        int ir = 0;
        for (int ix = 0; ix < this->rho_basis->nx; ++ix)
        {
            const Treal phase_x = kq_pack->kv_ptr->kvec_d[ikq].x * ix / this->rho_basis->nx;
            for (int iy = 0; iy < this->rho_basis->ny; ++iy)
            {
                const Treal phase_xy = phase_x + kq_pack->kv_ptr->kvec_d[ikq].y * iy / this->rho_basis->ny;
                for (int iz = this->rho_basis->startz_current; iz < this->rho_basis->startz_current + this->rho_basis->nplane; ++iz)
                {
                    const Treal phase_xyz = phase_xy + kq_pack->kv_ptr->kvec_d[ikq].z * iz / this->rho_basis->nz;
                    const T exp_tmp = std::exp(phase_xyz * this->two_pi_i);
                    this->phi[iw][ir] = porter[ir] * exp_tmp;
                    ++ir;
                }
            }
        }
    }
    ModuleBase::timer::tick("Exx_Lip", "phi_cal");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::psi_cal()
{
    ModuleBase::TITLE("Exx_Lip", "psi_cal");
    ModuleBase::timer::tick("Exx_Lip", "psi_cal");
    if (PARAM.inp.init_chg == "atomic")
    {
        std::vector<T> porter (this->wfc_basis->nrxx);
        for (int iq = 0; iq < this->q_pack->kv_ptr->get_nks(); ++iq)
        {
            for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
            {
                this->wfc_basis->recip2real(&(this->q_pack->kspw_psi_ptr->operator()(iq, ib, 0)), porter.data(), iq);

                int ir = 0;
                for (int ix = 0; ix < this->rho_basis->nx; ++ix)
                {
                    const Treal phase_x = this->q_pack->kv_ptr->kvec_d[iq].x * ix / this->rho_basis->nx;
                    for (int iy = 0; iy < this->rho_basis->ny; ++iy)
                    {
                        const Treal phase_xy = phase_x + this->q_pack->kv_ptr->kvec_d[iq].y * iy / this->rho_basis->ny;
                        for (int iz = this->rho_basis->startz_current; iz < this->rho_basis->startz_current + this->rho_basis->nplane; ++iz)
                        {
                            const Treal phase_xyz = phase_xy + this->q_pack->kv_ptr->kvec_d[iq].z * iz / this->rho_basis->nz;
                            const T exp_tmp = std::exp(phase_xyz * this->two_pi_i);
                            this->psi[iq][ib][ir] = porter[ir] * exp_tmp;
                            ++ir;
                        }
                    }
                }
            }
        }
    }
    else if (PARAM.inp.init_chg == "file")
    {
        for (int iq = 0; iq < this->q_pack->kv_ptr->get_nks(); ++iq)
        {
            this->phi_cal(this->q_pack, iq);
            for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
            {
                ModuleBase::GlobalFunc::ZEROS(this->psi[iq][ib].data(), this->rho_basis->nrxx);
                for (int iw = 0; iw < PARAM.globalv.nlocal; ++iw)
                {
                    for (int ir = 0; ir < this->rho_basis->nrxx; ++ir)
                    {
                        this->psi[iq][ib][ir] += (*this->q_pack->hvec_array)(iq, ib, iw) * this->phi[iw][ir];
                    }
                }
            }
        }
    }
    ///////////////////////////////////////////////////////////////////////////
    //			!!! k_point parallel incompleted. need to loop iq in other pool
    ///////////////////////////////////////////////////////////////////////////
    ModuleBase::timer::tick("Exx_Lip", "psi_cal");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::judge_singularity(const int ik)
{
    ModuleBase::timer::tick("Exx_Lip", "judge_singularity");
    if (PARAM.inp.init_chg=="atomic")
    {
        this->iq_vecik = ik;
    }
    else if(PARAM.inp.init_chg=="file")
    {
        Treal min_q_minus_k(std::numeric_limits<Treal>::max());
        for( int iq=0; iq<this->q_pack->kv_ptr->get_nks(); ++iq)
        {
            const Treal q_minus_k((this->q_pack->kv_ptr->kvec_c[iq] - this->k_pack->kv_ptr->kvec_c[ik]).norm2());
            if(q_minus_k < min_q_minus_k)
            {
                min_q_minus_k = q_minus_k;
                this->iq_vecik = iq;
            }
        }
    }
    ModuleBase::timer::tick("Exx_Lip", "judge_singularity");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::qkg2_exp(const int ik, const int iq)
{
    ModuleBase::timer::tick("Exx_Lip", "qkg2_exp");
    for( int ig=0; ig<this->rho_basis->npw; ++ig)
    {
        const Treal qkg2 = ((this->q_pack->kv_ptr->kvec_c[iq] - this->k_pack->kv_ptr->kvec_c[ik] + this->rho_basis->gcar[ig]) * (ModuleBase::TWO_PI / this->ucell_ptr->lat0)).norm2();
        if (Conv_Coulomb_Pot_K::Ccp_Type::Ccp == info.ccp_type || Conv_Coulomb_Pot_K::Ccp_Type::Hf == info.ccp_type)
        {
            if (std::abs(qkg2) < 1e-10)
                { this->recip_qkg2[ig] = 0.0; }												// 0 to ignore bb/qkg2 when qkg2==0
            else
                { this->recip_qkg2[ig] = 1.0 / qkg2; }
            this->sum2_factor += this->recip_qkg2[ig] * std::exp(-info.lambda * qkg2);
            this->recip_qkg2[ig] = sqrt(this->recip_qkg2[ig]);
        }
        else if (Conv_Coulomb_Pot_K::Ccp_Type::Hse == info.ccp_type)
        {
            if (std::abs(qkg2) < 1e-10)
                { this->recip_qkg2[ig] = 1.0 / (2 * info.hse_omega); }
            else
                { this->recip_qkg2[ig] = sqrt((1 - std::exp(-qkg2 / (4 * info.hse_omega * info.hse_omega))) / qkg2); }
        }
    }
    ModuleBase::timer::tick("Exx_Lip", "qkg2_exp");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::b_cal(const int ik, const int iq, const int ib)
{
    ModuleBase::timer::tick("Exx_Lip", "b_cal");
    const ModuleBase::Vector3<double> q_minus_k = this->q_pack->kv_ptr->kvec_d[iq] - this->k_pack->kv_ptr->kvec_d[ik];
    std::vector<T > mul_tmp(this->rho_basis->nrxx);
    for( size_t ir=0,ix=0; ix<this->rho_basis->nx; ++ix)
    {
        const Treal phase_x = q_minus_k.x * ix / this->rho_basis->nx;
        for( size_t iy=0; iy<this->rho_basis->ny; ++iy)
        {
            const Treal phase_xy = phase_x + q_minus_k.y * iy / this->rho_basis->ny;
            for( size_t iz=this->rho_basis->startz_current; iz<this->rho_basis->startz_current+this->rho_basis->nplane; ++iz)
            {
                const Treal phase_xyz = phase_xy + q_minus_k.z * iz / this->rho_basis->nz;
                mul_tmp[ir] = std::exp(-phase_xyz * this->two_pi_i);
                mul_tmp[ir] *= this->psi[iq][ib][ir];
                ++ir;
            }
        }
    }

    std::vector<T> porter (this->rho_basis->nrxx);
    for(size_t iw=0; iw< PARAM.globalv.nlocal; ++iw)
    {
        auto& phi_w = this->phi[iw];
        for( size_t ir=0; ir<this->rho_basis->nrxx; ++ir)
        {
            porter[ir] = conj(phi_w[ir]) * mul_tmp[ir] ;
            // porter[ir] = phi_w[ir] * psi_q_b[ir] *exp_tmp[ir] ;
        }
        T* const b_w = &this->b[iw * this->rho_basis->npw];
        this->rho_basis->real2recip( porter.data(), b_w);
        if (Conv_Coulomb_Pot_K::Ccp_Type::Ccp == info.ccp_type || Conv_Coulomb_Pot_K::Ccp_Type::Hf == info.ccp_type) {
            if ((iq == this->iq_vecik) && (gzero_rank_in_pool == GlobalV::RANK_IN_POOL)) {							/// need to check while use k_point parallel
                this->b0[iw] = b_w[this->rho_basis->ig_gge0];
        } }

        for (size_t ig = 0; ig < this->rho_basis->npw; ++ig)
            { b_w[ig] *= this->recip_qkg2[ig]; }
    }
    ModuleBase::timer::tick("Exx_Lip", "b_cal");
}

template <typename T, typename Device>
void  Exx_Lip<T, Device>::sum3_cal(const int iq, const int ib)
{
    ModuleBase::timer::tick("Exx_Lip", "sum3_cal");
    if (gzero_rank_in_pool == GlobalV::RANK_IN_POOL) {
        for (int iw_l = 0; iw_l < PARAM.globalv.nlocal; ++iw_l) {
            for (int iw_r = 0; iw_r < PARAM.globalv.nlocal; ++iw_r) {
                this->sum3[iw_l][iw_r] += this->b0[iw_l] * conj(this->b0[iw_r]) * (Treal)this->q_pack->wf_wg(iq, ib);
    } } }
    ModuleBase::timer::tick("Exx_Lip", "sum3_cal");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::b_sum(const int iq, const int ib)			// Peize Lin change 2019-04-14
{
    ModuleBase::timer::tick("Exx_Lip", "b_sum");
    // this->sum1[iw_l,iw_r] += \sum_{ig} this->b[iw_l,ig] * conj(this->b[iw_r,ig]) * this->q_pack->wf_wg(iq,ib)
    LapackConnector::herk(
        'U','N',
        PARAM.globalv.nlocal, this->rho_basis->npw,
        (Treal)this->q_pack->wf_wg(iq, ib), this->b.data(), this->rho_basis->npw,
        1.0, this->sum1.data(), PARAM.globalv.nlocal);
    // cblas_zherk( CblasRowMajor, CblasUpper, CblasNoTrans,
    // 			PARAM.globalv.nlocal, this->rho_basis->npw,
    // 			this->q_pack->wf_wg(iq,ib), static_cast<void*>(this->b), this->rho_basis->npw,
    // 			1.0, static_cast<void*>(this->sum1), PARAM.globalv.nlocal);
    ModuleBase::timer::tick("Exx_Lip", "b_sum");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::sum_all(const int ik)
{
    ModuleBase::timer::tick("Exx_Lip", "sum_all");
    Treal sum2_factor_g = 0.0;
    const Treal fourpi_div_omega = 4 * (Treal)(ModuleBase::PI / this->ucell_ptr->omega);
    const Treal spin_fac = 2.0;
  #ifdef __MPI
    if (Conv_Coulomb_Pot_K::Ccp_Type::Ccp == info.ccp_type || Conv_Coulomb_Pot_K::Ccp_Type::Hf == info.ccp_type)
        { MPI_Reduce(&this->sum2_factor, &sum2_factor_g, 1, MPI_DOUBLE, MPI_SUM, gzero_rank_in_pool, POOL_WORLD); }
  #endif
    for (size_t iw_l = 1; iw_l < PARAM.globalv.nlocal; ++iw_l) {
        for (size_t iw_r = 0; iw_r < iw_l; ++iw_r) {
            this->sum1[iw_l * PARAM.globalv.nlocal + iw_r] = conj(this->sum1[iw_r * PARAM.globalv.nlocal + iw_l]);		// Peize Lin add conj 2019-04-14
    } }

    for (int iw_l = 0; iw_l < PARAM.globalv.nlocal; ++iw_l)
    {
        for( int iw_r=0; iw_r<PARAM.globalv.nlocal; ++iw_r)
        {
            this->exx_matrix[ik][iw_l][iw_r] = -fourpi_div_omega * this->sum1[iw_l * PARAM.globalv.nlocal + iw_r] * spin_fac;
            if (Conv_Coulomb_Pot_K::Ccp_Type::Ccp == info.ccp_type || Conv_Coulomb_Pot_K::Ccp_Type::Hf == info.ccp_type)
            {
                if (gzero_rank_in_pool == GlobalV::RANK_IN_POOL)
                {
                    this->exx_matrix[ik][iw_l][iw_r] += spin_fac * (fourpi_div_omega * this->sum3[iw_l][iw_r] * sum2_factor_g);
                    this->exx_matrix[ik][iw_l][iw_r] += spin_fac * (-1 / (Treal)sqrt(info.lambda * ModuleBase::PI) * (Treal)(this->q_pack->kv_ptr->get_nks() / PARAM.inp.nspin) * this->sum3[iw_l][iw_r]);
                }
            }
        }
    }
    ModuleBase::timer::tick("Exx_Lip", "sum_all");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::exx_energy_cal()
{
    ModuleBase::TITLE("Exx_Lip","exx_energy_cal");
    ModuleBase::timer::tick("Exx_Lip", "exx_energy_cal");

    Treal exx_energy_tmp = 0.0;

    for( int ik=0; ik<this->k_pack->kv_ptr->get_nks(); ++ik) {
        for( int iw_l=0; iw_l<PARAM.globalv.nlocal; ++iw_l) {
            for( int iw_r=0; iw_r<PARAM.globalv.nlocal; ++iw_r) {
                for( int ib=0; ib<PARAM.inp.nbands; ++ib) {
                    exx_energy_tmp += (this->exx_matrix[ik][iw_l][iw_r] * conj((*this->k_pack->hvec_array)(ik, ib, iw_l)) * (*this->k_pack->hvec_array)(ik, ib, iw_r)).real() * this->k_pack->wf_wg(ik, ib);
    } } } }
  #ifdef __MPI
    MPI_Allreduce( &exx_energy_tmp, &this->exx_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);				// !!! k_point parallel incompleted. different pools have different kv.set_nks(>) deadlock
  #endif
    this->exx_energy *= (PARAM.inp.nspin==1) ? 2 : 1;
    this->exx_energy /= 2;										// ETOT = E_band - 1/2 E_exx
    ModuleBase::timer::tick("Exx_Lip", "exx_energy_cal");
}

template <typename T, typename Device>
void Exx_Lip<T, Device>::write_q_pack() const
{
    ModuleBase::timer::tick("Exx_Lip", "write_q_pack");

    if (PARAM.inp.out_chg[0] == 0)
        { return; }

    if (!GlobalV::RANK_IN_POOL)
    {
        const std::string exx_q_pack = "exx_q_pack/";

        const std::string command_mkdir = "test -d " + PARAM.globalv.global_out_dir + exx_q_pack + " || mkdir " + PARAM.globalv.global_out_dir + exx_q_pack;
        assert( system(command_mkdir.c_str()) == 0);

        const std::string command_kpoint = "test -f " + PARAM.globalv.global_out_dir + exx_q_pack + PARAM.inp.kpoint_file + " || cp " + PARAM.inp.kpoint_file + " " + PARAM.globalv.global_out_dir + exx_q_pack + PARAM.inp.kpoint_file;
        assert( system(command_kpoint.c_str()) == 0);

        std::stringstream ss_wf_wg;
        ss_wf_wg << PARAM.globalv.global_out_dir << exx_q_pack << "wf_wg_" << GlobalV::MY_POOL;
        std::ofstream ofs_wf_wg(ss_wf_wg.str().c_str());
        for( int iq = 0; iq < this->q_pack->kv_ptr->get_nks(); ++iq)
        {
            for( int ib=0; ib<PARAM.inp.nbands; ++ib)
            {
                ofs_wf_wg<<this->q_pack->wf_wg(iq,ib)<<"\t";
            }
            ofs_wf_wg<<std::endl;
        }
        ofs_wf_wg.close();

        std::stringstream ss_hvec;
        ss_hvec	<< PARAM.globalv.global_out_dir << exx_q_pack << "hvec_" << GlobalV::MY_POOL;
        std::ofstream ofs_hvec(ss_hvec.str().c_str());
        for( int iq=0; iq<this->q_pack->kv_ptr->get_nks(); ++iq)
        {
            for( int iw=0; iw<PARAM.globalv.nlocal; ++iw)
            {
                for( int ib=0; ib<PARAM.inp.nbands; ++ib)
                {
                    ofs_hvec << (*this->q_pack->hvec_array)(iq, ib, iw).real() << " " << (*this->q_pack->hvec_array)(iq, ib, iw).imag() << " ";
                }
                ofs_hvec<<std::endl;
            }
        }
        ofs_hvec.close();
    }
    ModuleBase::timer::tick("Exx_Lip", "write_q_pack");
}

/*
void Exx_Lip::read_q_pack(const ModuleSymmetry::Symmetry& symm,
                          const ModulePW::PW_Basis_K* this->wfc_basis,
                          const Structure_Factor& sf)
{
    const std::string exx_q_pack = "exx_q_pack/";
    this->q_pack = new k_package();
    this->q_pack->kv_ptr = new K_Vectors();
    const std::string exx_kpoint_card = PARAM.globalv.global_out_dir + exx_q_pack + PARAM.inp.kpoint_file;
    this->q_pack->kv_ptr->set( symm, exx_kpoint_card, PARAM.inp.nspin, this->ucell_ptr->G, this->ucell_ptr->latvec, GlobalV::ofs_running );
    this->q_pack->wf_ptr = new wavefunc();
    this->q_pack->wf_ptr->allocate(this->q_pack->kv_ptr->get_nkstot(),
                             this->q_pack->kv_ptr->get_nks(),
                             this->q_pack->kv_ptr->ngk.data(),
                             this->wfc_basis->npwk_max); // mohan update 2021-02-25
    //	this->q_pack->wf_ptr->init(this->q_pack->kv_ptr->get_nks(),this->q_pack->kv_ptr,this->ucell_ptr,old_pwptr,&ppcell,&GlobalC::ORB,&hm,&Pkpoints);
    this->q_pack->wf_ptr->table_local.create(GlobalC::ucell.ntype, GlobalC::ucell.nmax_total, PARAM.globalv.nqx);
    // this->q_pack->wf_ptr->table_local.create(this->q_pack->wf_ptr->this->ucell_ptr->ntype, this->q_pack->wf_ptr->this->ucell_ptr->nmax_total, PARAM.globalv.nqx);
  #ifdef __LCAO
    Wavefunc_in_pw::make_table_q(GlobalC::ORB.orbital_file, this->q_pack->wf_ptr->table_local);
    // Wavefunc_in_pw::make_table_q(this->q_pack->wf_ptr->ORB_ptr->orbital_file, this->q_pack->wf_ptr->table_local, this->q_pack->wf_ptr);
    for(int iq=0; iq<this->q_pack->kv_ptr->get_nks(); ++iq)
    {
        Wavefunc_in_pw::produce_local_basis_in_pw(iq,
                                                  this->wfc_basis,
                                                  sf,
                                                  this->q_pack->wf_ptr->wanf2[iq],
                                                  this->q_pack->wf_ptr->table_local);
        //		Wavefunc_in_pw::produce_local_basis_in_pw(iq, this->q_pack->wf_ptr->wanf2[iq], this->q_pack->wf_ptr->table_local,
        // this->q_pack->wf_ptr);
    }
  #endif
    this->q_pack->wf_wg.create(this->q_pack->kv_ptr->get_nks(),PARAM.inp.nbands);
    if(!GlobalV::RANK_IN_POOL)
    {
        std::stringstream ss_wf_wg;
        ss_wf_wg << PARAM.globalv.global_out_dir << exx_q_pack << "wf_wg_" << GlobalV::MY_POOL;
        std::ifstream ifs_wf_wg(ss_wf_wg.str().c_str());
        for( int iq = 0; iq < this->q_pack->kv_ptr->get_nks(); ++iq)
        {
            for( int ib=0; ib<PARAM.inp.nbands; ++ib)
            {
                ifs_wf_wg>>this->q_pack->wf_wg(iq,ib);
            }
        }
        ifs_wf_wg.close();
    }
    #ifdef __MPI
    MPI_Bcast( this->q_pack->wf_wg.c, this->q_pack->kv_ptr->get_nks()*PARAM.inp.nbands, MPI_DOUBLE, 0, POOL_WORLD);
    #endif
    this->q_pack->hvec_array = new ModuleBase::ComplexMatrix [this->q_pack->kv_ptr->get_nks()];
    for( int iq=0; iq<this->q_pack->kv_ptr->get_nks(); ++iq)
    {
        this->q_pack->hvec_array[iq].create(PARAM.globalv.nlocal,PARAM.inp.nbands);
    }
    if(!GlobalV::RANK_IN_POOL)
    {
        std::stringstream ss_hvec;
        ss_hvec	<< PARAM.globalv.global_out_dir << exx_q_pack << "hvec_" << GlobalV::MY_POOL;
        std::ifstream ifs_hvec(ss_hvec.str().c_str());
        for( int iq=0; iq<this->q_pack->kv_ptr->get_nks(); ++iq)
        {
            for( int iw=0; iw<PARAM.globalv.nlocal; ++iw)
            {
                for( int ib=0; ib<PARAM.inp.nbands; ++ib)
                {
                    double a,this->b;
                    ifs_hvec>>a>>this->b;
                    this->q_pack->hvec_array[iq](iw,ib) = {a,this->b};
                }
            }
        }
        ifs_hvec.close();
    }
    #ifdef __MPI
    for( int iq=0; iq<this->q_pack->kv_ptr->get_nks(); ++iq)
    {
        MPI_Bcast( this->q_pack->hvec_array[iq].c, PARAM.globalv.nlocal*PARAM.inp.nbands, MPI_DOUBLE_COMPLEX, 0, POOL_WORLD);
    }
    #endif
    return;
}
*/

#endif