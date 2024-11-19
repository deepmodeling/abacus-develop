#include <iostream>

#include "module_base/matrix.h"
#include "module_base/name_angular.h"
#include "module_base/scalapack_connector.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"
//#include "module_hamilt_pw/hamilt_pwdft/onsite_projector.h"
#include "spin_constrain.h"
#include "module_parameter/parameter.h"
#ifdef __LCAO
#include "module_elecstate/elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/operator_lcao/dspin_lcao.h"

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mi_lcao(const int& step, bool print)
{
    ModuleBase::TITLE("module_deltaspin", "cal_mi_lcao");
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "cal_mi_lcao");
    // calculate MW from lambda in real space projection method
    this->zero_Mi();
    const hamilt::HContainer<double>* dmr
        = static_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM()->get_DMR_pointer(1);
    std::vector<double> moments;
    if(PARAM.inp.nspin==2)
    {
        static_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM()->switch_dmr(2);
        moments = static_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, double>>*>(this->p_operator)->cal_moment(dmr, this->get_constrain());
        static_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM()->switch_dmr(0);
        for(int iat=0;iat<this->Mi_.size();iat++)
        {
            this->Mi_[iat].x = 0.0;
            this->Mi_[iat].y = 0.0;
            this->Mi_[iat].z = moments[iat];
        }
    }
    else if(PARAM.inp.nspin==4)
    {
        moments = static_cast<hamilt::DeltaSpin<hamilt::OperatorLCAO<std::complex<double>, std::complex<double>>>*>(this->p_operator)->cal_moment(dmr, this->get_constrain());
        for(int iat=0;iat<this->Mi_.size();iat++)
        {
            this->Mi_[iat].x = moments[iat*3];
            this->Mi_[iat].y = moments[iat*3+1];
            this->Mi_[iat].z = moments[iat*3+2];
        }
    }

    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "cal_mi_lcao");
}

#endif

/*template <>
void spinconstrain::SpinConstrain<std::complex<double>>::cal_mi_pw()
{
    ModuleBase::TITLE("module_deltaspin", "cal_mi_pw");
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "cal_mi_pw");

    this->zero_Mi();
    if(PARAM.inp.device == "cpu")
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_CPU>::get_instance();
        // loop over k-points to calculate Mi of \sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
        std::complex<double>* psi_pointer = nullptr;
        psi::Psi<std::complex<double>, base_device::DEVICE_CPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_CPU>*>(this->psi);
        const int nbands = psi_t->get_nbands();
        const int nks = psi_t->get_nk();
        const int npol = psi_t->npol;
        for(int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik); // tabulate for each atom at each k-point
            // std::cout << __FILE__ << ":" << __LINE__ << " nbands = " << nbands << std::endl;
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer);
            const std::complex<double>* becp = onsite_p->get_h_becp();
            // becp(nbands*npol , nkb)
            // mag = wg * \sum_{nh}becp * becp
            int nkb = onsite_p->get_tot_nproj();
            for(int ib = 0;ib<nbands;ib++)
            {
                const double weight = this->pelec->wg(ik, ib);
                int begin_ih = 0;
                for(int iat = 0; iat < this->Mi_.size(); iat++)
                {
                    std::complex<double> occ[4] = {ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO};
                    const int nh = onsite_p->get_nh(iat);
                    for(int ih = 0; ih < nh; ih++)
                    {
                        const int index = ib*2*nkb + begin_ih + ih;
                        occ[0] += conj(becp[index]) * becp[index];
                        occ[1] += conj(becp[index]) * becp[index + nkb];
                        occ[2] += conj(becp[index + nkb]) * becp[index];
                        occ[3] += conj(becp[index + nkb]) * becp[index + nkb];
                    }
                    // occ has been reduced and calculate mag
                    this->Mi_[iat].z += weight * (occ[0] - occ[3]).real();
                    this->Mi_[iat].x += weight * (occ[1] + occ[2]).real();
                    this->Mi_[iat].y += weight * (occ[1] - occ[2]).imag();
                    begin_ih += nh;
                }
            }
        }
    }
#if ((defined __CUDA) || (defined __ROCM))
    else
    {
        auto* onsite_p = projectors::OnsiteProjector<double, base_device::DEVICE_GPU>::get_instance();
        // loop over k-points to calculate Mi of \sum_{k,i,l,m}<Psi_{k,i}|alpha_{l,m}><alpha_{l,m}|Psi_{k,i}>
        std::complex<double>* psi_pointer = nullptr;
        psi::Psi<std::complex<double>, base_device::DEVICE_GPU>* psi_t = static_cast<psi::Psi<std::complex<double>, base_device::DEVICE_GPU>*>(this->psi);
        const int nbands = psi_t->get_nbands();
        const int nks = psi_t->get_nk();
        const int npol = psi_t->npol;
        for(int ik = 0; ik < nks; ik++)
        {
            psi_t->fix_k(ik);
            psi_pointer = psi_t->get_pointer();
            onsite_p->tabulate_atomic(ik); // tabulate for each atom at each k-point
            // std::cout << __FILE__ << ":" << __LINE__ << " nbands = " << nbands << std::endl;
            onsite_p->overlap_proj_psi(nbands * npol, psi_pointer);
            const std::complex<double>* becp = onsite_p->get_h_becp();
            // becp(nbands*npol , nkb)
            // mag = wg * \sum_{nh}becp * becp
            int nkb = onsite_p->get_size_becp() / nbands / npol;
            for(int ib = 0;ib<nbands;ib++)
            {
                const double weight = this->pelec->wg(ik, ib);
                int begin_ih = 0;
                for(int iat = 0; iat < this->Mi_.size(); iat++)
                {
                    std::complex<double> occ[4] = {ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO};
                    const int nh = onsite_p->get_nh(iat);
                    for(int ih = 0; ih < nh; ih++)
                    {
                        const int index = ib*2*nkb + begin_ih + ih;
                        occ[0] += conj(becp[index]) * becp[index];
                        occ[1] += conj(becp[index]) * becp[index + nkb];
                        occ[2] += conj(becp[index + nkb]) * becp[index];
                        occ[3] += conj(becp[index + nkb]) * becp[index + nkb];
                    }
                    // occ has been reduced and calculate mag
                    this->Mi_[iat].z += weight * (occ[0] - occ[3]).real();
                    this->Mi_[iat].x += weight * (occ[1] + occ[2]).real();
                    this->Mi_[iat].y += weight * (occ[1] - occ[2]).imag();
                    begin_ih += nh;
                }
            }
        }
    }
#endif
    // reduce mag from all k-pools
    Parallel_Reduce::reduce_double_allpool(PARAM.inp.kpar, GlobalV::NPROC_IN_POOL, &(this->Mi_[0][0]), 3 * this->Mi_.size());
    
    ModuleBase::timer::tick("spinconstrain::SpinConstrain", "cal_mi_pw");
}*/

template <>
void spinconstrain::SpinConstrain<std::complex<double>>::set_operator(
    hamilt::Operator<std::complex<double>>* op_in)
{
    this->p_operator = op_in;
}

template <>
void spinconstrain::SpinConstrain<double>::set_operator(
    hamilt::Operator<double>* op_in)
{
    this->p_operator = op_in;
}