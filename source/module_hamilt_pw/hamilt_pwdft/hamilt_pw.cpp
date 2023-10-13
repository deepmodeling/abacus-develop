#include "hamilt_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include "operator_pw/veff_pw.h"
#include "operator_pw/ekinetic_pw.h"
#include "operator_pw/meta_pw.h"
#include "operator_pw/nonlocal_pw.h"

namespace hamilt
{

template<typename T, typename Device>
HamiltPW<T, Device>::HamiltPW(elecstate::Potential* pot_in, ModulePW::PW_Basis_K* wfc_basis, K_Vectors* pkv)
{
    this->classname = "HamiltPW";
    const auto tpiba2 = static_cast<Real>(GlobalC::ucell.tpiba2);
    const auto tpiba = static_cast<Real>(GlobalC::ucell.tpiba);
    const int* isk = pkv->isk.data();
    const Real* gk2 = wfc_basis->get_gk2_data<Real>();

    if (GlobalV::T_IN_H)
    {
        // Operator<double>* ekinetic = new Ekinetic<OperatorLCAO<double>>
        Operator<T, Device>* ekinetic
            = new Ekinetic<OperatorPW<T, Device>>(tpiba2, gk2, wfc_basis->nks, wfc_basis->npwk_max);
        if(this->ops == nullptr)
        {
            this->ops = ekinetic;
        }
        else
        {
            this->ops->add(ekinetic);
        }
    }
    if (GlobalV::VL_IN_H)
    {
        std::vector<std::string> pot_register_in;
        if (GlobalV::VION_IN_H)
        {
            pot_register_in.push_back("local");
        }
        if (GlobalV::VH_IN_H)
        {
            pot_register_in.push_back("hartree");
        }
        //no variable can choose xc, maybe it is necessary
        pot_register_in.push_back("xc");
        if (GlobalV::imp_sol)
        {
            pot_register_in.push_back("surchem");
        }
        if (GlobalV::EFIELD_FLAG)
        {
            pot_register_in.push_back("efield");
        }
        if (GlobalV::GATE_FLAG)
        {
            pot_register_in.push_back("gatefield");
        }
        //only Potential is not empty, Veff and Meta are available
        if(pot_register_in.size()>0)
        {
            //register Potential by gathered operator
            pot_in->pot_register(pot_register_in);
            Operator<T, Device>* veff
                = new Veff<OperatorPW<T, Device>>(isk,
                                                       pot_in->get_v_effective_data<Real>(),
                                                       pot_in->get_effective_v().nr,
                                                       pot_in->get_effective_v().nc,
                                                       wfc_basis);
            if(this->ops == nullptr)
            {
                this->ops = veff;
            }
            else
            {
                this->ops->add(veff);
            }
            Operator<T, Device>* meta
                = new Meta<OperatorPW<T, Device>>(tpiba,
                                                       isk,
                                                       pot_in->get_vofk_effective_data<Real>(),
                                                       pot_in->get_effective_vofk().nr,
                                                       pot_in->get_effective_vofk().nc,
                                                       wfc_basis);
            this->ops->add(meta);
        }
    }
    if (GlobalV::VNL_IN_H)
    {
        Operator<T, Device>* nonlocal
            = new Nonlocal<OperatorPW<T, Device>>(isk, &GlobalC::ppcell, &GlobalC::ucell, wfc_basis);
        if(this->ops == nullptr)
        {
            this->ops = nonlocal;
        }
        else
        {
            this->ops->add(nonlocal);
        }
    }
    return;
}

template<typename T, typename Device>
HamiltPW<T, Device>::~HamiltPW()
{
    if(this->ops!= nullptr)
    {
        delete this->ops;
    }
}

template<typename T, typename Device>
void HamiltPW<T, Device>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltPW","updateHk");
    this->ops->init(ik);
    ModuleBase::TITLE("HamiltPW","updateHk");
}

template<typename T, typename Device>
template<typename T_in, typename Device_in>
HamiltPW<T, Device>::HamiltPW(const HamiltPW<T_in, Device_in> *hamilt)
{
    this->classname = hamilt->classname;
    OperatorPW<std::complex<T_in>, Device_in> * node =
            reinterpret_cast<OperatorPW<std::complex<T_in>, Device_in> *>(hamilt->ops);

    while(node != nullptr) {
        if (node->classname == "Ekinetic") {
            Operator<T, Device>* ekinetic =
                    new Ekinetic<OperatorPW<T, Device>>(
                            reinterpret_cast<const Ekinetic<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = ekinetic;
            }
            else {
                this->ops->add(ekinetic);
            }
            // this->ops = reinterpret_cast<Operator<T, Device>*>(node);
        }
        else if (node->classname == "Nonlocal") {
            Operator<T, Device>* nonlocal =
                    new Nonlocal<OperatorPW<T, Device>>(
                            reinterpret_cast<const Nonlocal<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = nonlocal;
            }
            else {
                this->ops->add(nonlocal);
            }
        }
        else if (node->classname == "Veff") {
            Operator<T, Device>* veff =
                    new Veff<OperatorPW<T, Device>>(
                            reinterpret_cast<const Veff<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = veff;
            }
            else {
                this->ops->add(veff);
            }
        }
        else if (node->classname == "Meta") {
            Operator<T, Device>* meta =
                    new Meta<OperatorPW<T, Device>>(
                            reinterpret_cast<const Meta<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = meta;
            }
            else {
                this->ops->add(meta);
            }
        }
        else {
            ModuleBase::WARNING_QUIT("HamiltPW", "Unrecognized Operator type!");
        }
        node = reinterpret_cast<OperatorPW<std::complex<T_in>, Device_in> *>(node->next_op);
    }
}

template class HamiltPW<std::complex<float>, psi::DEVICE_CPU>;
template class HamiltPW<std::complex<double>, psi::DEVICE_CPU>;
// template HamiltPW<std::complex<double>, psi::DEVICE_CPU>::HamiltPW(const HamiltPW<std::complex<double>, psi::DEVICE_CPU> *hamilt);
#if ((defined __CUDA) || (defined __ROCM))
template class HamiltPW<std::complex<float>, psi::DEVICE_GPU>;
template class HamiltPW<std::complex<double>, psi::DEVICE_GPU>;
// template HamiltPW<std::complex<double>, psi::DEVICE_GPU>::HamiltPW(const HamiltPW<std::complex<double>, psi::DEVICE_GPU> *hamilt);
#endif

} // namespace hamilt